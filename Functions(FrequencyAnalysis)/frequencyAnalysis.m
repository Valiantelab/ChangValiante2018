function [epileptiformEvent, frequencyContentAnalysis] = frequencyAnalysis (InputGUI, PathName, FileName)
%% Stage 2: Re-process the File (after manual changes) for further analysis
% Author: Michael Chang
% Run this file after the detection algorithm to analyze the results and do
% additional analysis to the detected events. This creates the time vector,
% LFP time series, LED if there is light, and filters the data using a
% bandpass filter (1-50 Hz) and a low pass filter (@68 Hz). It will also
% report the dominant frequency content of a time series.

%User Input variables
userInput = str2double(InputGUI(1:3)); %convert inputs into numbers

    % Load .mat file
    matFileName = fullfile(PathName,FileName);  %location of the .matfile    
    myVars = {'artifactSpikes', 'details', 'samplingInterval', 'SLE', 'spikes', 'x', 'metadata'};   %variables of interest from .mat file
    load(sprintf('%s', matFileName), myVars{:}) %load .mat file into workspace
    
    % Auto-load .abf file; to load 'x' into workspace (can't save if it's too large)
    abfFileName = sprintf('%s.abf', matFileName(1:end-4));  %Locate the .abf file from the same folder as the .mat file     
    x=abfload(abfFileName); %load 'x' from .abf into workspace
    
% %Load excel file - with human adjusted onset/offset times
%     [FileName,PathName] = uigetfile ('*.xls','pick .xls file to load excel sheet', inputdir);%Choose file    
%     excel_filename = fullfile(PathName, FileName);
    
    %Auto-load excel file with onset/offset times into workspace
    excel_filename = sprintf('%s(organized).xls', matFileName(1:end-4));    %name of .xls file in the same folder as .mat file
    
    excelSheet = xlsread(excel_filename, 3); %Read the 3rd sheets in excel file
    
    events = excelSheet(2:end,:);   %events selected after editing the excel sheet by humans

f = waitbar(0,'Loading Data: Please wait while data is prepared...','Name', 'Stage 2 Analysis: in progress');

%Create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%Seperate signals from .abf files
LFP = x(:,1);   %original LFP signal
if size (x,2) > 1
    LED = x(:,size(x,2));   %light pulse signal, assumed to be from the last channel
    onsetDelayLED = 0.13;  %seconds
    offsetDelayLED = 1.5;  %seconds
    lightpulse = LED > 1;
else
    LED =[];
    onsetDelay = [];
end

% Find Light pulse
if ~isempty(LED)       
    [P] = pulse_seq(LED);   %determine location of light pulses     

    %Find range of time when light pulse has potential to trigger an event,
    clear lightTriggeredRange lightTriggeredZone
    for i = 1:numel(P.range(:,1))
        lightTriggeredOnsetRange = (P.range(i,1):P.range(i,1)+(onsetDelayLED*frequency));
        lightTriggeredOnsetZone{i} = lightTriggeredOnsetRange; 
        lightTriggeredOffsetRange = (P.range(i,1):P.range(i,1)+(offsetDelayLED*frequency));
        lightTriggeredOffsetZone{i} = lightTriggeredOffsetRange; 
    end

    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredOnsetZones = cat(2, lightTriggeredOnsetZone{:});  %2 is vertcat
    lightTriggeredOffsetZones = cat(2, lightTriggeredOffsetZone{:});  %2 is vertcat
end


%% Filter Bank
waitbar(0.02, f, 'Please wait while Processing Data...');

%Center the LFP data
LFP_centered = LFP - LFP(1);

%Bandpass butter filter [1 - 100 Hz]
[b,a] = butter(2, ([1 100]/(frequency/2)), 'bandpass');
LFP_filtered = filtfilt (b,a,LFP);             %Filtered signal

%Absolute value of the filtered data
AbsLFP_Filtered = abs(LFP_filtered);            %1st derived signal

%Derivative of the filtered data (absolute value)
DiffLFP_Filtered = abs(diff(LFP_filtered));     %2nd derived signal

%Power of the filtered data (feature for classification)
powerFeature = (LFP_filtered).^2;                     %3rd derived signal
avgPowerFeature = mean(powerFeature);   %for use as the intensity ratio threshold, later


%% Part 2 - Feature Extraction: Duration, Intensity, and Peak-to-Peak Amplitude
waitbar(0.1, f, 'Extracting Features from detected Epileptiform Events');
intensityPerMinute = cell(size(events,1),1);
totalPower = zeros(size(events,1),1);
for i = 1:size(events,1)
    %Convert times (s) into points in data
    onsetTime = int64(events(i,1)*frequency);
    offsetTime = int64(events(i,2)*frequency);
    
    %Incase crawler function sets offsetTime prior to onsetTime
    if onsetTime > offsetTime   %onset comes after offset, means error has occured
        events(i,2) = events(i,1) + 1; %add 1 second to onset to find arbitary offset
        offsetTime = int64(events(i,2)*frequency);
    end
    
    %make epileptiform event vector
    eventVector = int64(onsetTime:offsetTime);  %SLE Vector

    %Split the event vector into (1 min) windows
    window_size = 1;  %seconds; you can change window size as desired
    sleDuration = round(numel(eventVector)/frequency);    %rounded to whole number; Note: sometimes the SLE crawler can drop the duration of the event to <1 s
    if sleDuration == 0
        sleDuration = 1;    %need to fix this so you don't analyze event vectors shorter than 1 s
        fprintf(2,'\nWarning! You detected a epileptiform that is shorter than 1 sec, this is an issue with your SLECrawler.m algorithm.')   %testing if I still need this if statement, I think I've fixed teh algorithm so no events <1 s have their features extracted
    end

    %Calculate the intensity (per sec) for epileptiform events; intensity (per sec) is average power
    clear spikeRateMinute intensity    
    for j = 1:sleDuration
        startWindow = onsetTime+((window_size*frequency)*(j-1));          
        endWindow = onsetTime+((window_size*frequency)*j);
        %Calculate the intensity per minute for epileptiform events
        if numel(powerFeature) > endWindow
            PowerPerMinute = sum(powerFeature (startWindow:endWindow));
        else
            PowerPerMinute = sum(powerFeature (startWindow:numel(powerFeature)));
        end
        intensity(j,1) = startWindow; %time windows starts
        intensity(j,2) = PowerPerMinute;   %Total power within the (minute) window
    end

    intensityPerMinute{i} = intensity;    %store the intensity per minute of each SLE for analysis later
    
    %Calculate average intensity of epileptiform event
    totalPower(i) = sum(powerFeature(eventVector)); %This is the energy of the event; totalPower is a misnomer
    events (i,5) = totalPower(i) /sleDuration;  %This is the power of the event; you previously considered it as the average power/minute

    %Calculate peak-to-peak amplitude of epileptiform event
    eventVectorLFP = LFP_centered(eventVector);
    p2pAmplitude = max(eventVectorLFP) - min (eventVectorLFP);
    events (i,6) = p2pAmplitude;

    %Calculate the event's duration 
    events (i,3) = events (i,2) - events (i,1);
    
    %m calculation
%     test=WP_MultipleRegression(eventVectorLFP', 10000);

end

%Calculate dominant frequency
waitbar(0.35, f, 'Calculating Frequency Content of Epileptiform Events');
pptx_filename = sprintf('%s(frequencyContent).pptx', matFileName(1:end-4));    %name of .pptx file in the same folder as .mat file
[epileptiformEvent, interictal] = dominantFrequency(spikes, events, SLE, artifactSpikes, samplingInterval, x, pptx_filename, userInput(1), userInput(2), userInput(3));    %userInput(3) is figure 0=no 1=yes

%Determine if event is light-triggered, delay to onset
waitbar(0.55, f, 'Determine if Epileptiform Events are light-triggered');
if ~isempty(LED)
    %Classify which SLEs were light triggered | if the spike onset is after light pulse
    for i=1:size(events,1) 
        %use the "ismember" function 
        events(i,8)=ismember (int64(events(i,1)*frequency), lightTriggeredOnsetZones);    
    end
    
    %Back-up algorithm to detect which SLEs were light triggered | if the peak of spike is after light pulse; in case the crawler function detects the onset prior to light pulse
    for i=1:size(events,1) 
        if events(i,8) == 1
            continue
        else            
        %use the "ismember" function 
        events(i,8)=ismember (int64(events(i,28)*frequency), lightTriggeredOnsetZones);    
        end
    end
         
    %Find the delay between ictal event onset (peak of the first spike) and all light pulses
    calcIctalOnsetDelays(numel(P.range(:,1)),1) = 0; %Preallocate 
    
    %Calculate total number of light pulses
    totalLightPulses = numel(P.range(:,1)); %Permanent code line
    
    %Calculate delay between ictal event onset (peak of the first spike) and preceeding light pulse 
    for i=1:size(events,1) 
        clear calcIctalOnsetDelays
        %Calculate the delay between peak of first spike in ictal event and all light pulses                           
        for j = 1:numel(P.range(:,1))
            calcIctalOnsetDelays(j) = events(i,28) - P.range(j,1)/frequency;
        end
        %Determine location of the preceding light pulse
        if any(calcIctalOnsetDelays>0)   %In case a spontaneous ictal event occurs before the first light pulse, there will be no preceding light pulse
            [ictalOnsetDelay, index_ictalOnsetDelay] = min(calcIctalOnsetDelays(calcIctalOnsetDelays>0));  %Find the smallest positive value
        else
            ictalOnsetDelay = 0;    %0 = spontaneous event, no preceding light pulse; didn't want to put Nan in case it threw off my code.
            index_ictalOnsetDelay = 1;  %Arbitarily use the first light pulse to push the code forwards
        end
        
        %Calculate the interstimulus interval
        if index_ictalOnsetDelay < totalLightPulses %In the event the last stimulus triggers an event, there won't be a subseqent stimulus afterwardsto to subtract
            interstimulusInterval = (P.range(index_ictalOnsetDelay+1,1) - P.range(index_ictalOnsetDelay,1))/frequency;
        else
            interstimulusInterval = (numel(LFP) - P.range(index_ictalOnsetDelay,1))/frequency;  %consider the interstimulus interval to end when the recording ends
        end
        %Store the ictal onset delay and interstimulus interval
        events (i, 26) = ictalOnsetDelay;
        events (i, 27) = interstimulusInterval;                
    end
end

%% Stage 3: Analyze the file
% Author: Michael Chang
% Run this file after Stage 2 to analyze the results statistically

waitbar(0.75, f, 'Statistical Analysis of Epileptiform Events');

%organize groups index
indexControl = events(:,4)==1;
indexTest = events(:,4)==2;
indexPosttest = events(:,4)==3;

%% Duration
waitbar(0.80, f, 'Statistical Analysis of Epileptiform Events: Duration');
%Organize groups 
feature = 3;    %column #, i.e., duration is the 3rd column
durationControl = events(events(:,4)==1,feature);
durationTest = events(events(:,4)==2,feature);
durationPosttest = events(events(:,4)==3,feature);

%duration Matrix
durationMatrix(1:numel(durationControl),1) = durationControl;
durationMatrix(1:numel(durationTest),2) = durationTest;
durationMatrix(1:numel(durationPosttest),3) = durationPosttest;
durationMatrix(durationMatrix==0) = NaN;

%Analysis
%Control Condition
[resultsDuration(1,:)] = stage3Analysis (durationControl, 'nonparametric', 'no figure');
%Test Condition
resultsDuration(2,:) = stage3Analysis (durationTest, 'nonparametric', 'no figure');
%Posttest Conditions
if numel(durationPosttest) >2
resultsDuration(3,:) = stage3Analysis (durationPosttest, 'nonparametric', 'no figure');   
end

%Comparisons
[h,p,D] = kstest2(durationControl, durationTest); %2-sample KS Test
p_value_KS_duration = p;    %KS Test's p value
d_KS_duration = D;  %KS Test's D statistic
cliffs_d_duration = CliffDelta(durationControl, durationTest);   %Cliff's d

% %Mann Whitney U Test (aka Wilcoxin Ranked Sum Test)
% [p,h,stats] = ranksum(durationControl, durationTest)

% %Independent sample Student's T-test
% [h,p,ci,stats] = ttest2(durationControl, durationTest);

%one-way ANOVA
[p,tbl,stats] = anova1(durationMatrix);
tbl_1ANOVA_duration = tbl;

% title('Box Plot of ictal event duration from different time periods')
% xlabel ('Time Period')
% ylabel ('Duration of ictal events (s)')

%multiple comparisons, Tukey-Kramer Method
c = multcompare(stats);
c_duration = c;

%Kruskal Wallis
% p = kruskalwallis(durationMatrix);
% title('Boxplot: duration of ictal events from different treatment groups')
% xlabel ('Treatment Group')
% ylabel ('Duration (s)')


%% Intensity
waitbar(0.85, f, 'Statistical Analysis of Epileptiform Events: Intensity');
%Organize groups
feature = 5;
intensityControl = events(events(:,4)==1,feature);
intensityTest = events(events(:,4)==2,feature);
intensityPosttest = events(events(:,4)==3,feature);

%intensity Matrix
intensityMatrix(1:numel(intensityControl),1) = intensityControl;
intensityMatrix(1:numel(intensityTest),2) = intensityTest;
intensityMatrix(1:numel(intensityPosttest),3) = intensityPosttest;
intensityMatrix(intensityMatrix==0) = NaN;

%Analysis
%Control Condition
resultsIntensity(1,:) = stage3Analysis (intensityControl, 'nonparametric', 'nofigure');
%Test Condition
resultsIntensity(2,:) = stage3Analysis (intensityTest, 'nonparametric', 'nofigure');
%Posttest Conditions
if numel(intensityPosttest)>2
resultsIntensity(3,:) = stage3Analysis (intensityPosttest, 'nonparametric', 'nofigure');
end

%Analysis, comparison
%2-sample KS Test
[h,p,D] = kstest2 (intensityControl, intensityTest);    %KS Test's p value
p_value_KS_intensity = p; %KS Test's D statistic
d_KS_intensity = D; %KS Test's D statistic
cliffs_d_intensity = CliffDelta(intensityControl,intensityTest); %Cliff's D

% %Independent sample Student's T-test
% [h,p,ci,stats] = ttest2(intensityControl, intensityTest);

%one-way ANOVA
[p,tbl,stats] = anova1(intensityMatrix);
tbl_1ANOVA_intensity = tbl;

% title('Box Plot of ictal event intensity from different time periods')
% xlabel ('Treatment Condition')
% ylabel ('intensity of ictal events (mV^2/s)')

%Multiple Comparisons, Tukey-Kramer Method
c = multcompare(stats);
c_intensity = c;

% %Kruskal Wallis
% p = kruskalwallis(intensityMatrix);
% title('Boxplot intensity of ictal events from different treatment groups')
% xlabel ('Treatment Group')
% ylabel ('intensity (e-5)')

%% Circular Variance and Plots of Ictal Event with photosimulation    
waitbar(0.83, f, 'Statistical Analysis of Epileptiform Events: Circ Stat');
%Calculate Theta
for i = 1:numel(events(:,1))
    events(i, 29) = events(i, 26)/events(i, 27) * (2*pi);
end
%Organize Group
feature = 29;
thetaControl = events(events(:,4)==1,feature);
thetaTest = events(events(:,4)==2,feature);
thetaPosttest = events(events(:,4)==3,feature);

% thetaControl=SLE(controlStart<SLE(:,1) & SLE(:,1)<controlEnd,feature);
% thetaTest=SLE(testStart<SLE(:,1) & SLE(:,1)<testEnd,feature);
% thetaPosttest=SLE(posttestStart<SLE(:,1) & SLE(:,1)<posttestEnd,feature);

% Analysis, light correlation?
resultsTheta(1,1)=circ_vtest(thetaControl,0);
resultsTheta(2,1)=circ_vtest(thetaTest,0);

%Figures for Visual Analysis
if userInput(3) == 1
    FigE=figure;
    set(gcf,'Name','Control', 'NumberTitle', 'off');
    circ_plot(thetaControl,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Control Condition, p = %.3f', resultsTheta(1,1)));

    FigF=figure;
    set(gcf,'Name','Test','NumberTitle', 'off');
    circ_plot(thetaTest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Test Condition, p = %.3f', resultsTheta(2,1)));
end

if numel(thetaPosttest)>2
    resultsTheta(3,1)=circ_vtest(thetaPosttest,0);
    %Figures for visual analysis
    if userInput(3) == 1
        FigG=figure;
        set(gcf,'Name','Post-Test','NumberTitle', 'off');
        circ_plot(thetaPosttest,'hist',[],50,false,true,'linewidth',2,'color','r');
        title (sprintf('Post-Test Condition, p = %.3f',resultsTheta(3,1)));
    end    
end

%% Dominant Frequency Content; the frequency that carries more power (PSD) w.r.t.
%https://dsp.stackexchange.com/questions/40180/the-exact-definition-of-dominant-frequency
waitbar(0.85, f, 'Statistical Analysis of Epileptiform Events: Dominant Frequency');

 %Preallocate
frequencyContentAnalysis = zeros(size(epileptiformEvent,1),9);

%Re-calculate window size (for frequency analysis); previous window size
%was for analyzing 'intensity'; terrible mistake naming the different
%vaiable with the same name

% windowOverlap = epileptiformEvent{i,2}(1);    %the first window size represents the overlap
% windowSize = windowOverlap*2;
windowSize = userInput(1);

indexEvents = find(events(:,3) > windowSize);
% for i = 1:size(epileptiformEvent,1)
for i = indexEvents'

    %Place time next to the dominant frequency for each epileptiform event
    epileptiformEvent{i,3}(:,2) = epileptiformEvent{i,2};
    
    %Calculate Tonic Phase
    tonicPhase = 5; %Hz
    dominantFreq = epileptiformEvent{i,3}(:,1:2);
    indexTonic = dominantFreq(:,1) > tonicPhase;   %locate the period of time when ictal event has tonic phase
    epileptiformEvent{i,3}(:,3) = indexTonic;  %store Boolean index          
    
    %locate start of Tonic phase | Contingous segments above threshold    
    for j = 1: numel (indexTonic) %slide along the SLE; 
        if j+1 > numel(indexTonic)  %If you scan through the entire SLE and don't find a tonic phase, classify event as a IIE                 
            j = find(indexTonic(:),1,'first');  %Take the 1st sec where frequency is 'high' as the onset if back-to-back high frequency are not found
                if isempty(j) %There is no Tonic Phase
                    j = 1;  %honory position, just to push the function through 
%                     classification = 0;       
                else if j == numel(indexTonic)
                        j = j-1;
                    else
                        j = j;
                    end
                end                
            startTonicPhase(1,1) = dominantFreq(j,2);    %time (s)
            startTonicPhase(1,2) = j;    %store the index                
%                 while ~and(indexTonic(j) == 0, indexTonic(j+1) == 0) & dominantFreq(j,1) ~= 0; %Locate the offset time as the first point as either two continueous segments with low frequency or zero frequency, which ever comes first
                while ~and(indexTonic(j) == 0, indexTonic(j+1) == 0) & j+1 < numel(indexTonic); %Locate the offset time as the first point as either two continueous segments with low frequency or zero frequency, which ever comes first
                    j = j+1;    %keep sliding along the SLE until the statement above is false.
                    if j+1 > numel(indexTonic)  %If you slide all the way to the end and still can't find tonic phase offset, 
                        j = numel(indexTonic)+1;  %take the last point as the offset of the tonic phase - this means there is no clonic phase; add 1 because note you will remove it in line 33
                        classification = 2;   %1 = SLE;   2 = Tonic-only     
                        break
                    end                                    
                end
            if j > 1    
                endTonic(1,1) = dominantFreq(j-1,2); %take the point before the frequency drops to be the end of the tonic-phase
            else
                j = 2;  %in the event, j = 1 to push the algorithm through (line 400), make j = 2 so you can find an arbitary endTonic point
                endTonic(1,1) = dominantFreq(j-1,2); %This is an arbitary point to push the algorithm though.
            end
                endTonic(1,2) = j-1;  %store the index  
            classification = 0; %There is no tonic-clonic ictal event, just forced the algorithm through to locate the potential tonic phase
        else                        
            if indexTonic(j) > 0 && indexTonic(j+1) > 0 %If you locate two continuous segments with high frequency, mark the first segment as start of tonic phase                        
                startTonicPhase(1,1) = dominantFreq(j,2);  %store the onset time
                startTonicPhase(1,2) = j;    %store the index
                classification = 1;   %1 = tonic-clonic SLE;   2 = tonic-only
                while ~and(indexTonic(j) == 0, indexTonic(j+1) == 0) & dominantFreq(j,1) ~= 0; %Locate the offset time as the first point as either two continueous segments with low frequency or zero frequency, which ever comes first
                    j = j+1;    %keep sliding along the SLE until the statement above is false.
                    if j+1 > numel(indexTonic)  %If you slide all the way to the end and still can't find tonic phase offset, 
                        j = numel(indexTonic)+1;  %take the last point as the offset of the tonic phase - this means there is no clonic phase; add 1 because note you will remove it in line 33
                        classification = 2;   %1 = SLE;   2 = Tonic-only     
                        break
                    end                                    
                end            
                endTonic(1,1) = dominantFreq(j-1,2); %take the point before the frequency drops to be the end of the tonic-phase
                endTonic(1,2) = j-1;  %store the index   
                break
                %This is to confirm the tonic phase is >3s
%                 if (endTonic(1)-startTonicPhase(1)) < 0.1  %There needs to be at least 3 seconds of tonic phase
%                     continue
%                 else
%                     break
%                 end                
            end
        end        
    end
    
    %Calculate average Frequency during tonic phase
    meanTonicFreq = mean(epileptiformEvent{i,3}(startTonicPhase(2):endTonic(2),1));

    %Calculate average Frequency during clonic phase
    meanClonicFreq = mean(epileptiformEvent{i,3}(endTonic(2)+1:end,1));   
    
    %Ictal Event Classification 
    frequencyContentAnalysis(i,1) = classification; %0 = no ictal; 1 = tonic-clonic ictal; 2 = tonic-only ictal
    %Ictal Event Onset Time
    frequencyContentAnalysis(i,2) = startTonicPhase(1)-windowSize; %remove the padding added on to ictal event onset
    %Ictal Evenet Offset Time
    frequencyContentAnalysis(i,3) = endTonic(1)-windowSize; %remove the padding added on to ictal event onset
    %Ictal Evenet Duration
    frequencyContentAnalysis(i,4) = endTonic(1) - startTonicPhase(1);        
    %Percentage into SLE tonic phase starts
    frequencyContentAnalysis(i,5) = frequencyContentAnalysis(i,2)/ events(i,3);        
    
    [~, maxIndex] = max(dominantFreq(:,1)); %Find index max frequency occurs
    %Max Frequency
    frequencyContentAnalysis(i,6) = dominantFreq(maxIndex,1);
    %Time the max Frequency occurs
    frequencyContentAnalysis(i,7) = dominantFreq(maxIndex,2) - windowSize;  

    %Mean Frequency during Tonic Phase
    frequencyContentAnalysis(i,8) = meanTonicFreq;

    %Mean Frequency during Clonic Phase
    frequencyContentAnalysis(i,9) = meanClonicFreq;
        
end

%organize groups
dominantFreqControl = frequencyContentAnalysis (indexControl, :);
dominantFreqTest = frequencyContentAnalysis (indexTest, :);
dominantFreqPosttest = frequencyContentAnalysis (indexPosttest, :);

%Calculate median values 
if sum(indexControl)>1
    medianFreqControl = median(dominantFreqControl);
else
    medianFreqControl = dominantFreqControl;
end

if sum(indexTest)>1    
    medianFreqTest = median(dominantFreqTest);
else
    medianFreqTest = (dominantFreqTest);
end

if sum(indexPosttest)>1
    medianFreqPosttest = median(dominantFreqPosttest);
else
    medianFreqPosttest = (dominantFreqPosttest);
end

%Concatenate for plotting later into excel sheets
medianDominantFreq = vertcat(medianFreqControl,medianFreqTest,medianFreqPosttest);

%Max Dominant Frequency Matrix
dominantFreqMatrix(1:size(dominantFreqControl,1),1) = dominantFreqControl(:,6);
dominantFreqMatrix(1:size(dominantFreqTest,1),2) = dominantFreqTest(:,6);
dominantFreqMatrix(1:size(dominantFreqPosttest,1),3) = dominantFreqPosttest(:,6);
dominantFreqMatrix(dominantFreqMatrix==0) = NaN;

%Analysis
%Control Condition
resultsDominantFreq(1,:) = stage3Analysis (dominantFreqControl(:,6), 'nonparametric', 'no figure');
%Test Condition
resultsDominantFreq(2,:) = stage3Analysis (dominantFreqTest(:,6), 'nonparametric', 'no figure');
%Posttest Conditions
if numel(dominantFreqPosttest(:,6)) >2
resultsDominantFreq(3,:) = stage3Analysis (dominantFreqPosttest(:,6), 'nonparametric', 'no figure');   
end

%Comparisons
[h,p,D] = kstest2(dominantFreqControl(:,6), dominantFreqTest(:,6)); %2-sample KS Test
p_value_KS_dominantFreq = p;    %KS Test's p value
d_KS_dominantFreq = D;  %KS Test's D statistic
cliffs_d_dominantFreq = CliffDelta(dominantFreqControl(:,6), dominantFreqTest(:,6));   %Cliff's d

%one-way ANOVA
[p,tbl,stats] = anova1(dominantFreqMatrix);
tbl_1ANOVA_dominantFreq = tbl;

%multiple comparisons, Tukey-Kramer Method
c = multcompare(stats);
c_dominantFreq = c;

%% Combine all the results
result = horzcat(resultsDuration(:,1:3),resultsIntensity(:,1:3),resultsTheta, resultsDominantFreq(:,1:3));  %only the first three columns

%ictal events # in each group
n(1,1)=numel(thetaControl);
n(2,1)=numel(thetaTest);
n(3,1)=numel(thetaPosttest);

%% Write results to .xls 
waitbar(0.90, f, 'Write Results to excel sheets');

% Customize Analysis 
excelFileName = InputGUI{5}; %excel sheet output is written
sheetName = FileName(1:8);

%Label treatment groups
treatmentGroups(1,1) = {'Control'};
if InputGUI{4} == ""
    treatmentGroups(2,1) = {'Test'};
else
    treatmentGroups(2,1) = {InputGUI{4}};
end
treatmentGroups(3,1) = {'PostTest'};
% treatmentGroups = [1:3]';


%set subtitle
A = 'Treatment Group';
B = 'Duration (s), median';
C = 'Duration (s), IQR';
D = 'AD test, normality';
E = 'Power (mV^2/s), median';
F = 'Power (mV^2/s), IQR';
G = 'AD test, normality';
H = 'Light-triggered';
I = 'Dominant Frequency, median';

II = 'Dominant Frequency, IQR';
JJ = 'AD test, normality';
KK = 'n';

J = 'KS Test, 1 vs 2';
K = 'one-way ANOVA, duration';
M = 'Multiple Comparison (Tukey-Kramer method), duration';
N = 'one-way ANOVA, power';
O = 'Multiple Comparison (Tukey-Kramer method), power';
NN = 'one-way ANOVA, dominant frequency';
OO = 'Multiple Comparison (Tukey-Kramer method), dominant frequency';

P = 'Group';
Q = 'p-value';

R = 'p-value';
S = 'KS D stat';
T = 'Cliffs D';

%% Write General Results for Statistical Analysis
    subtitle1 = {A, B, C, D, E, F, G, H, I, II, JJ, KK};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A1');
    xlswrite(sprintf('%s',excelFileName),treatmentGroups,sprintf('%s',sheetName),'A2');
    xlswrite(sprintf('%s',excelFileName),result,sprintf('%s',sheetName),'B2');    
    xlswrite(sprintf('%s',excelFileName),n,sprintf('%s',sheetName),'L2');
%Write KS Test results
    subtitle1 = {J};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A6');    
    xlswrite(sprintf('%s',excelFileName),p_value_KS_duration,sprintf('%s',sheetName),'B6'); %KS Test p-value
    xlswrite(sprintf('%s',excelFileName),d_KS_duration,sprintf('%s',sheetName),'C6');   %KS Test D Value
    xlswrite(sprintf('%s',excelFileName),cliffs_d_duration,sprintf('%s',sheetName),'D6');    %Cliff's D        
    xlswrite(sprintf('%s',excelFileName),p_value_KS_intensity,sprintf('%s',sheetName),'E6'); %KS Test p-value
    xlswrite(sprintf('%s',excelFileName),d_KS_intensity,sprintf('%s',sheetName),'F6'); %KS Test D value
    xlswrite(sprintf('%s',excelFileName),cliffs_d_intensity, sprintf('%s',sheetName),'G6'); %Cliff's D
    xlswrite(sprintf('%s',excelFileName),p_value_KS_dominantFreq,sprintf('%s',sheetName),'I6'); %KS Test p-value
    xlswrite(sprintf('%s',excelFileName),d_KS_dominantFreq,sprintf('%s',sheetName),'J6'); %KS Test D value
    xlswrite(sprintf('%s',excelFileName),cliffs_d_dominantFreq, sprintf('%s',sheetName),'K6'); %Cliff's D
    subtitle2 = {R,S,T};
    xlswrite(sprintf('%s',excelFileName),subtitle2,sprintf('%s',sheetName),'B5');
    xlswrite(sprintf('%s',excelFileName),subtitle2,sprintf('%s',sheetName),'E5');
    xlswrite(sprintf('%s',excelFileName),subtitle2,sprintf('%s',sheetName),'I5');

%Write one-way ANOVA results, duration
    subtitle1 = {K};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A8');
    xlswrite(sprintf('%s',excelFileName),tbl_1ANOVA_duration,sprintf('%s',sheetName),'A9');
%Write multiple comparison (Tukey-Kramer), duration
    subtitle1 = {M};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A14');
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'A15'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'B15'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{Q},sprintf('%s',sheetName),'F15'); %p-value subtitle
    xlswrite(sprintf('%s',excelFileName),c_duration,sprintf('%s',sheetName),'A16');
%Write one-way ANOVA results, intensity
    subtitle1 = {N};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A20');
    xlswrite(sprintf('%s',excelFileName),tbl_1ANOVA_intensity,sprintf('%s',sheetName),'A21');
%Write multiple comparison (Tukey-Kramer), duration
    subtitle1 = {O};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A26');
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'A27'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'B27'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{Q},sprintf('%s',sheetName),'F27'); %p-value subtitle
    xlswrite(sprintf('%s',excelFileName),c_intensity,sprintf('%s',sheetName),'A28');    
%Write one-way ANOVA results, dominant frequency
    subtitle1 = {NN};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A32');
    xlswrite(sprintf('%s',excelFileName),tbl_1ANOVA_dominantFreq,sprintf('%s',sheetName),'A33');
%Write multiple comparison (Tukey-Kramer), dominant frequency
    subtitle1 = {OO};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A38');
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'A39'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'B39'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{Q},sprintf('%s',sheetName),'F39'); %p-value subtitle
    xlswrite(sprintf('%s',excelFileName),c_dominantFreq,sprintf('%s',sheetName),'A40');

%% Write Results of frequency content analysis 
    FileName = excel_filename;  %rename excel filename
    
    subtitle1 = {'Onset (s)', 'Offset (s)', 'Duration (s)', A};
    xlswrite(sprintf('%s',FileName),subtitle1,'Freq Content Analysis','A1');    
    xlswrite(sprintf('%s',FileName),events(:,1:4),'Freq Content Analysis','A2');   
    subtitle1 = {'Classification', 'Onset (s), Tonic Phase', 'Offset (s), Tonic Phase',	'Duration (s), Tonic Phase', '% into SLE', 'Max Frequency (Hz)', 'Time of Max Frequency (s)', 'Tonic Freq (Hz), mean', 'Clonic Freq (Hz), mean'};
    xlswrite(sprintf('%s',FileName),subtitle1,'Freq Content Analysis','E1');
    xlswrite(sprintf('%s',FileName),frequencyContentAnalysis,'Freq Content Analysis','E2');   

%Write average dominant frequency content of epileptiform events
    subtitle1 = {A};    %Label Treatment Group
    xlswrite(sprintf('%s',FileName),subtitle1,'Median Freq Content','A1');    
    subtitle1 = {'Classification', 'Onset (s), Tonic Phase', 'Offset (s), Tonic Phase',	'Duration (s), Tonic Phase', '% into SLE', 'Max Frequency (Hz)', 'Time of Max Frequency (s)', 'Tonic Freq (Hz), mean', 'Clonic Freq (Hz), mean', KK};
    xlswrite(sprintf('%s',FileName),subtitle1,'Median Freq Content','B1');   
    xlswrite(sprintf('%s',FileName),treatmentGroups,'Median Freq Content','A2');   
    xlswrite(sprintf('%s',FileName),medianDominantFreq,'Median Freq Content','B2');
    xlswrite(sprintf('%s',FileName),n,'Median Freq Content','K2');   
    
%% save final output

close all %Close all figures
waitbar(0.93, f, 'Saving workspace(stage 2) as .mat file');
%     save(sprintf('%s(stage2).mat', FileName(1:end-4)))  %Save Workspace
    
waitbar(1, f, 'Stage 2 Analysis: Complete')

close (f)

fprintf(1,'\nStage 2 Analysis Complete: %s\n', FileName)
end


