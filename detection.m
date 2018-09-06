%Program: Epileptiform Activity Detector 
%Author: Michael Chang (michael.chang@live.ca), Fred Chen and Liam Long; 
%Copyright (c) 2018, Valiante Lab
%Version 4.0

%% Clear All
close all
clear all
clc

%% GUI to set thresholds
%Settings, request for user input on threshold
titleInput = 'Specify Detection Thresholds';
prompt1 = 'Epileptiform Spike Threshold: average + (3.5 x Sigma)';
prompt2 = 'Artifact Threshold: average + (100 x Sigma) ';
prompt3 = 'Figure: Yes (1) or No (0)';
prompt4 = 'Stimulus channel (enter 0 if none):';
prompt5 = 'Plot all epileptiform events (Yes (1) or No (0):';
prompt6 = 'Classification Report';
prompt = {prompt1, prompt2, prompt3, prompt4, prompt5, prompt6};
dims = [1 70];
definput = {'3.5', '100', '0', '2', '1', '0'};
opts = 'on';
userInput = str2double(inputdlg(prompt,titleInput,dims,definput, opts));

%setting on distance between spikes, hard coded
distanceSpike = 0.15;  %distance between spikes (seconds)
distanceArtifact = 0.6; %distance between artifacts (seconds)
minSLEduration = 3; %seconds; %change to 5 s if any detection issues 


%% Load .abf and excel data
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', 'C:\Users\User\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)\1) abf files');%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k

%Label for titles
excelFileName = FileName(1:8);
uniqueTitle = '(epileptiformEvents)';

%% create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%% Seperate signals from .abf files
LFP = x(:,1);   %original LFP signal
if userInput(4)>0
    LED = x(:,userInput(4));   %light pulse signal, as defined by user's input via GUI
    onsetDelay = 0.13;  %seconds
else
    LED =[];
    onsetDelay = [];
end

%% Data Processing 
%Center the LFP data
LFP_normalized = LFP - LFP(1);                                      %centered signal at 0, y-axis

%Lowpass butter filter [2Hz]
fc = 2; % Cut off frequency
[b,a] = butter(2,fc/(frequency/2)); % Butterworth filter of order 2
LFP_normalizedLowPassFiltered = filtfilt(b,a,LFP_normalized); % Will be the filtered signal

%Bandpass butter filter [1 - 100 Hz]
[b,a] = butter(2, [[1 100]/(frequency/2)], 'bandpass');
LFP_normalizedFiltered = filtfilt (b,a,LFP_normalized);             %Filtered signal

%Absolute value of the filtered data
AbsLFP_normalizedFiltered = abs(LFP_normalizedFiltered);            %1st derived signal

%Derivative of the filtered data (absolute value)
DiffLFP_normalizedFiltered = abs(diff(LFP_normalizedFiltered));     %2nd derived signal

%Power of the filtered data (feature for classification)     
powerFeature = (LFP_normalizedFiltered).^2;                     %3rd derived signal

%% Detect potential events (epileptiform/artifacts) | Derivative Values
[epileptiformLocation, artifacts, locs_spike_1st] = detectEvents (DiffLFP_normalizedFiltered, frequency);

%remove potential events
for i = 1:size(epileptiformLocation,1)
AbsLFP_normalizedFiltered (epileptiformLocation (i,1):epileptiformLocation (i,2)) = [-1];
end

%remove artifacts
for i = 1:size(artifacts,1)
AbsLFP_normalizedFiltered (artifacts(i,1):artifacts(i,2)) = [-1];
end

%Isolate baseline recording
AbsLFP_normalizedFiltered (AbsLFP_normalizedFiltered == -1) = [];
AbsLFP_normalizedFilteredBaseline = AbsLFP_normalizedFiltered; %Rename

%Characterize baseline features from absolute value of the filtered data 
avgBaseline = mean(AbsLFP_normalizedFilteredBaseline); %Average
sigmaBaseline = std(AbsLFP_normalizedFilteredBaseline); %Standard Deviation

%% Detect events (epileptiform/artifacts) | Absolute Values
%Recreate the Absolute filtered LFP (1st derived signal) vector
AbsLFP_normalizedFiltered = abs(LFP_normalizedFiltered); %the LFP analyzed

%Define thresholds for detection, using inputs from GUI
minPeakHeight = avgBaseline+(userInput(1)*sigmaBaseline);      %threshold for epileptiform spike detection
minPeakDistance = distanceSpike*frequency;                              %minimum distance spikes must be apart
minArtifactHeight = avgBaseline+(userInput(2)*sigmaBaseline);  %threshold for artifact spike detection
minArtifactDistance = distanceArtifact*frequency;                       %minimum distance artifact spikes must be apart

%Detect events
[epileptiformLocation, artifactLocation, locs_spike_2nd] = detectEvents (AbsLFP_normalizedFiltered, frequency, minPeakHeight, minPeakDistance, minArtifactHeight, minArtifactDistance);

%% Feature extraction - Troubleshooting the classifier


%Preallocate
epileptiformLocation(:, 4:8) = 0;

%% Stage 1: Initial Classifier (duration)
%putativeIIS = epileptiformLocation (epileptiformLocation (:,3)<(minSLEduration*frequency),:); 
%putativeSLE = epileptiformLocation (epileptiformLocation (:,3)>=(minSLEduration*frequency),:); 

indexIIS = (epileptiformLocation (:,3)<(minSLEduration*frequency));
epileptiformLocation (indexIIS,7) = 3;   %3 = IIS; 0 = unclassified.

indexEvents = (epileptiformLocation (:,3)>=(minSLEduration*frequency));
putativeSLE = epileptiformLocation(indexEvents,:); 

%% Stage 2: Feature Extraction (Duration, Spiking Frequency and Intensity)

    if userInput(5) == 1   
        %% Creating powerpoint slide
        isOpen  = exportToPPTX();
        if ~isempty(isOpen),
            % If PowerPoint already started, then close first and then open a new one
            exportToPPTX('close');
        end

        exportToPPTX('new','Dimensions',[12 6], ...
            'Title','Epileptiform Event Detector V4.0', ...
            'Author','Michael Chang', ...
            'Subject','Automatically generated PPTX file', ...
            'Comments','This file has been automatically generated by exportToPPTX');

        exportToPPTX('addslide');
        exportToPPTX('addtext', 'Troubleshooting: Epileptiform Events detected', 'Position',[2 1 8 2],...
                     'Horiz','center', 'Vert','middle', 'FontSize', 36);
        exportToPPTX('addtext', sprintf('File: %s', FileName), 'Position',[3 3 6 2],...
                     'Horiz','center', 'Vert','middle', 'FontSize', 20);
        exportToPPTX('addtext', 'By: Michael Chang and Christopher Lucasius', 'Position',[4 4 4 2],...
                     'Horiz','center', 'Vert','middle', 'FontSize', 20);     

        exportToPPTX('addslide');
        exportToPPTX('addtext', 'Legend', 'Position',[0 0 4 1],...
                     'Horiz','left', 'Vert','middle', 'FontSize', 24);
        exportToPPTX('addtext', 'Epileptiform spike is average + 6*SD of the baseline', 'Position',[0 1 6 1],...
                     'Horiz','left', 'Vert','middle', 'FontSize', 14);
        exportToPPTX('addtext', 'Artifacts are average + 100*SD', 'Position',[0 2 5 1],...
                     'Horiz','left', 'Vert','middle', 'FontSize', 14);
        exportToPPTX('addtext', 'SLE onset is the first peak in power (minimum 1/3 of the max amplitude spike)', 'Position',[0 3 5 1],...
                     'Horiz','left', 'Vert','middle', 'FontSize', 14);
        exportToPPTX('addtext', 'SLE offset is when power returns below baseline/2', 'Position',[0 4 5 1],...
                     'Horiz','left', 'Vert','middle', 'FontSize', 14);
        exportToPPTX('addtext', 'Note: The event have only been shifted alone the y-axis to start at position 0', 'Position',[0 5 5 1],...
                     'Horiz','left', 'Vert','middle', 'FontSize', 16);      
    end
    

for i = 1:size(putativeSLE,1)   
    %make SLE vector
    onsetTime = putativeSLE(i,1);
    offsetTime = putativeSLE(i,2);
    eventIndex = (onsetTime:offsetTime);  %SLE Vector  
        
    %Calculate the spiking rate for epileptiform events
    windowSize = 1;  %seconds      
    sleDuration = round(numel(eventIndex)/frequency);    %rounded to whole number
    clear spikeRateMinute
    for j = 1:sleDuration
        startWindow = onsetTime+((windowSize*frequency)*(j-1));
        EndWindow = onsetTime+((windowSize*frequency)*j);
        spikeRate = and(startWindow<=locs_spike_2nd, EndWindow >=locs_spike_2nd);
        spikeRateMinute(j,1) = startWindow; %time windows starts
        spikeRateMinute(j,2) = sum(spikeRate(:));   %number of spikes in the window
    end
    
    spikeFrequency{i} = spikeRateMinute;    %store the spike frequency of each SLE for plotting later
    
    %average spike rate of SLE
    putativeSLE (i,4) = mean(spikeRateMinute(:,2));
          
    %average intensity of SLE
    totalPower = sum(powerFeature(eventIndex));
    putativeSLE (i,5) = totalPower /sleDuration;    
    
    %peak-to-peak amplitude
    eventVectorLFP = LFP_normalized(eventIndex);
    p2pAmplitude = max(eventVectorLFP) - min (eventVectorLFP);
    putativeSLE (i,6) = p2pAmplitude;
                
    %% Optional: plot vectors for Troubleshooting            
    if userInput(5) == 1   
        
        %set variables
        data1 = LFP_normalized; %Time series to be plotted 
        lightpulse = LED > 1;

        %make background vector
        if (onsetTime >= 50001 && (offsetTime+50000)<numel(data1))
            backgroundVector = (onsetTime-50000:offsetTime+50000);   %Background Vector
        elseif (onsetTime < 50001)
            backgroundVector = (1:offsetTime+50000);
        elseif ((offsetTime+50000)>numel(data1))
            backgroundVector = (onsetTime-50000:numel(data1));
        end
    
        %Plot figures
        figHandle = figure;
        set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
        set(gcf,'Name', sprintf ('Putative SLE #%d', i)); %select the name you want
        set(gcf, 'Position', get(0, 'Screensize'));   

        plot (t(backgroundVector),data1(backgroundVector))
        hold on
        plot (t(eventIndex),data1(eventIndex))     %SLE
        plot (t(onsetTime), data1(onsetTime), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
        plot (t(offsetTime), data1(offsetTime), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
        indexSpikes = and(onsetTime<locs_spike_2nd, offsetTime>locs_spike_2nd); %Locate spikes between the onset and offset  
        plot (t(locs_spike_2nd(indexSpikes)), (data1(locs_spike_2nd(indexSpikes))), 'x') %plot spikes (artifact removed)
        plot (t(backgroundVector),(lightpulse(backgroundVector)-1)/5, 'b') %plot LED   
        title (sprintf('LFP Recording, SLE #%d | For Troubleshooting', i));
        ylabel ('mV');
        xlabel ('Time (sec)');   

        yyaxis right

        plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'color', 'b')
        ylabel ('spike rate/second (Hz)');
        set(gca,'fontsize',20)

        exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
        exportToPPTX('addpicture',figHandle);      
        close(figHandle)
    end        
end

    if userInput(5) == 1   
    % save and close the .PPTX
    newFile = exportToPPTX('saveandclose',sprintf('%s%s', excelFileName, uniqueTitle)); 
    end
    
    
%% Stage 3: Final Classifier (k-means clustering)
%perform k-means clustering on the three feature sets

%classify based on average frequency 
[indexFrequency, thresholdFrequency] = sleClassifier (putativeSLE(:,4))   
putativeSLE (:,9) = putativeSLE (:,4)>thresholdFrequency;
    if userInput(6) == 1  
    %plot figure
    featureSet = putativeSLE(:,4);
    index = indexFrequency;
    featureThreshold = thresholdFrequency;
    michaelThreshold = 1 %Hz 
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Feature Set: Spiking Rate (Hz)'); %select the name you want
    gscatter(featureSet , featureSet, index);    %plot scatter plot
    hold on
    %plot the algorithm detected threshold
    plot ([featureThreshold featureThreshold], ylim); 
    %plot Michael Chang's threshold values 
    plot ([michaelThreshold michaelThreshold], ylim);
    %Label
    title ('Unsupervised classication, using k-means clustering');
    ylabel ('Spiking Rate (Hz)');
    xlabel ('Spiking Rate (Hz)');   
    legend('SLE', 'IIE', 'Algo Threshold', 'Michaels Threshold')
    set(gca,'fontsize',12)
    end


%classify based on average intensity 
[indexIntensity, thresholdIntensity] = sleClassifier (putativeSLE(:,5))
putativeSLE (:,10) = indexIntensity;

    if userInput(6) == 1  
    %plot figure
    featureSet = putativeSLE(:,5);
    index = indexIntensity;
    featureThreshold = thresholdIntensity;
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Feature Set: Intensity (Power/Duration)'); %select the name you want
    gscatter(featureSet , featureSet, index);    %plot scatter plot
    hold on
    %plot the algorithm detected threshold
    plot ([featureThreshold featureThreshold], ylim); 
    %plot Michael Chang's threshold values 
    if mean(featureSet)>std(featureSet)
        michaelThreshold = mean(featureSet)-std(featureSet);
    else 
        michaelThreshold = mean(featureSet);
    end
    plot ([michaelThreshold michaelThreshold], ylim);
    %Label
    title ('Unsupervised classication, using k-means clustering');
    ylabel ('Average Intensity (Power/Duration)');
    xlabel ('Average Intensity (Power/Duration)');   
    legend('SLE', 'IIE', 'Algo Threshold', 'Michaels Threshold')
    set(gca,'fontsize',12)
    end


%classify based on average amplitude 
[indexAmplitude, thresholdAmplitude] = sleClassifier (putativeSLE(:,6))   

    if userInput(6) == 1  
    %plot figure
    featureSet = putativeSLE(:,6);
    index = indexAmplitude;
    featureThreshold = thresholdAmplitude;
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Feature Set: Peak-to-Peak Amplitude (mV)'); %select the name you want
    gscatter(featureSet , featureSet, index);    %plot scatter plot
    hold on
    %plot the algorithm detected threshold
    plot ([featureThreshold featureThreshold], ylim); 
    %plot Michael Chang's threshold values 
    michaelThreshold = 100*std(featureSet);
    plot ([michaelThreshold michaelThreshold], ylim);
    %Label
    title ('Unsupervised classication, using k-means clustering');
    ylabel ('Peak-to-Peak Amplitude (mV)');
    xlabel ('Peak-to-Peak Amplitude (mV)');   
    legend('SLE', 'IIE', 'Algo Threshold', 'Michaels Threshold')
    set(gca,'fontsize',12)
    end

%Classifier 
for i = 1: numel(putativeSLE(:,1))
    if indexFrequency(i) + indexIntensity(i) +  indexAmplitude(i) == 3
    putativeSLE (i,7) = 1   %1 = SLE; 2 = IIE; 3 = IIS; 0 = unclassified.
    else
    putativeSLE (i,7) = 2
    end
end

indexSLE = putativeSLE (:,7) == 1;  %classify which ones are SLEs


% %Trouble shooting, temporary array variables 
% SLE_final = putativeSLE((putativeSLE (:,7) == 1),:);
% indexSLE = putativeSLE (:,7) == 1;
% SLE_final(:,1:3)=SLE_final(:,1:3)/frequency;
% IIS = epileptiformLocation (indexIIS,:);

    if userInput(6) == 1  
    %3D scatter plot
    figure;
    scatter3(putativeSLE(:,4), putativeSLE(:,5), putativeSLE(:,6), 18, 'red', 'filled')
    hold on
    scatter3(putativeSLE(indexSLE ,4), putativeSLE(indexSLE ,5), putativeSLE(indexSLE ,6), 18, 'blue', 'filled')
    %Label
    title ('Classified Epileptiform Events');
    xlabel ('Average Spiking Rate (Hz)');   
    ylabel ('Average Intensity (Power/Duration)');
    zlabel ('Peak-to-Peak Amplitude (mV)');    
    legend ('IIE', 'SLE')
    end


%% SLE: Determine exact onset and offset times | Power Feature
%Scan Low-Pass Filtered Power signal for precise onset/offset times
SLE_final = SLECrawler(LFP_normalizedFiltered, putativeSLE(indexSLE,1:3), frequency, LED, onsetDelay, locs_spike_2nd);  %can also define if light triggered


%Store light-triggered events (s)
%triggeredEvents = SLE_final(SLE_final(:,4)>0, :);

%% Write to .xls
%set subtitle
A = 'Onset (s)';
B = 'Offset (s)';
C = 'Duration (s)';
D = 'Light-triggered (1 = yes)';
E = 'Average Spike Rate (Hz)';
F = 'Average Intensity (power/duration)';

%Sheet 0 = Details - To be completed at a later date with Liam's help.
% https://www.mathworks.com/help/matlab/matlab_prog/creating-a-map-object.html
% details{1,1} = 'FileName:';     details {1,2} = sprintf('%s', FileName);
% details{2,1} = 'LED:';     details {2,2} = sprintf('%s', FileName);
% 
% details {2,1} = 'LED:';         
% 'Sampling Frequency:'
% 'Multiple of Sigma for spike threshold:' 
% 'Epileptiform Spike Threshold:'
% 'Minimum distance between epileptiform spikes:'
% 'Multiple of Sigma for artifact threshold:' 
% 'Artifact threshold:'
% 'Minimum distance between artifacts:'
% 'Minimum seizure duration:' 
% 'Maximum onset delay for stimulus'
% 
%     subtitle0 = {details(:,1)};
%     xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),subtitle0,'Details','A1');
%     xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),artifacts/frequency,'Artifacts','A2');
    

%Sheet 1 = Artifacts   
if isempty(artifactLocation) == 0
    subtitle3 = {A, B, C};
    xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),subtitle3,'Artifacts','A1');
    xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),artifactLocation/frequency,'Artifacts','A2');
else
    display ('No artifacts were detected.');
end

%Sheet 2 = IIS
if isempty(IIS) == 0  
    subtitle2 = {A, B, C, D};
    xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),subtitle2,'IIS' ,'A1');
    xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),IIS,'IIS' ,'A2');
else
    display ('No IISs were detected.');
end
    
%Sheet 3 = SLE
if isempty(SLE_final) == 0   
    subtitle1 = {A, B, C, D, E, F};
    xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),subtitle1,'SLE' ,'A1');
    xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),SLE_final,'SLE' ,'A2');
else
    display ('No SLEs were detected.');
end

%% Optional: Plot Figures
if userInput(3) == 1   
    
%% Creating powerpoint slide
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end

exportToPPTX('new','Dimensions',[12 6], ...
    'Title','Epileptiform Detector V4.0', ...
    'Author','Michael Chang', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');

exportToPPTX('addslide');
exportToPPTX('addtext', 'SLE Events detected ', 'Position',[2 1 8 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 36);
exportToPPTX('addtext', sprintf('File: %s', FileName), 'Position',[3 3 6 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
exportToPPTX('addtext', 'By: Michael Chang and Christopher Lucasius', 'Position',[4 4 4 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);     
         
exportToPPTX('addslide');
exportToPPTX('addtext', 'Legend', 'Position',[0 0 4 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 24);
exportToPPTX('addtext', 'Epileptiform spike is average + 6*SD of the baseline', 'Position',[0 1 6 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'Artifacts are average + 100*SD', 'Position',[0 2 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'SLE onset is the first peak in power (minimum 1/3 of the max amplitude spike)', 'Position',[0 3 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'SLE offset is when power returns below baseline/2', 'Position',[0 4 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'Note: The event have only been shifted alone the y-axis to start at position 0', 'Position',[0 5 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 16);      

%% plot entire recording 
figHandle = figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name', sprintf ('Overview of %s', FileName)); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

lightpulse = LED > 1;

subplot (3,1,1)
reduce_plot (t, LFP_normalized, 'k');
hold on
reduce_plot (t, lightpulse - 2);

%plot artifacts (red), found in 2nd search
for i = 1:numel(artifactLocation(:,1)) 
    reduce_plot (t(artifactLocation(i,1):artifactLocation(i,2)), LFP_normalized(artifactLocation(i,1):artifactLocation(i,2)), 'r');
end

%plot onset markers
for i=1:numel(epileptiformTime(:,1))
reduce_plot ((onsetTimes(i)), (LFP_normalized(epileptiformLocation(i))), 'o');
end

%plot offset markers
for i=1:numel(epileptiformTime(:,2))
reduce_plot ((offsetTimes(i)), (LFP_normalized(epileptiformLocation(i,2))), 'x');
end

title (sprintf ('Overview of LFP (10000 points/s), %s', FileName));
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,2) 
reduce_plot (t, AbsLFP_normalizedFiltered, 'b');
hold on
reduce_plot (t, lightpulse - 1);

%plot spikes (artifact removed)
for i=1:size(locs_spike_2nd,1)
plot (t(locs_spike_2nd(i,1)), (DiffLFP_normalizedFiltered(locs_spike_2nd(i,1))), 'x')
end

title ('Overview of filtered LFP (bandpass: 1 to 100 Hz)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,3) 
reduce_plot (t(1:end-1), DiffLFP_normalizedFiltered, 'g');
hold on

%plot spikes 
for i=1:size(locs_spike_1st,1)
plot (t(locs_spike_1st(i,1)), (DiffLFP_normalizedFiltered(locs_spike_1st(i,1))), 'x')
end

title ('Peaks (o) in Derivative of filtered LFP');
ylabel ('Derivative (mV)');
xlabel ('Time (s)');

exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
exportToPPTX('addpicture',figHandle);      
close(figHandle)

%% Plotting out detected SLEs with context | To figure out how off you are
data1 = LFP_normalized; %Time series to be plotted 

    lightpulse = LED > 1;

for i = 1:size(SLE_final,1)
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('V4.0 SLE #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
   
    onsetTime = single(SLE_final(i,1)*10000);
    offsetTime = single(SLE_final(i,2)*10000);
    sleVector = (onsetTime:offsetTime);  %SLE Vector    
    
    if (onsetTime >= 50001 && (offsetTime+50000)<numel(data1))
        backgroundVector = (onsetTime-50000:offsetTime+50000);   %Background Vector
    elseif (onsetTime < 50001)
        backgroundVector = (1:offsetTime+50000);
    elseif ((offsetTime+50000)>numel(data1))
        backgroundVector = (onsetTime-50000:numel(data1));
    end
    
    normalizeLFP = (data1(backgroundVector(1)));
    normalizeLED = abs(min(data1(sleVector)-normalizeLFP));
    plot (t(backgroundVector),data1(backgroundVector)-normalizeLFP ) %background
    hold on
    plot (t(backgroundVector),(lightpulse(backgroundVector)/4)-normalizeLED, 'b') %plot LED   
    plot (t(sleVector),data1(sleVector)-normalizeLFP)     %SLE
    plot (t(onsetTime), data1(onsetTime)-normalizeLFP , 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
    plot (t(offsetTime), data1(offsetTime)-normalizeLFP , 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
    indexSpikes = and(onsetTime<locs_spike_2nd, offsetTime>locs_spike_2nd); %Locate spikes between the onset and offset  
    plot (t(locs_spike_2nd(indexSpikes)), (data1(locs_spike_2nd(indexSpikes))-normalizeLFP), 'x', 'color', 'green') %plot spikes (artifact removed)
    
       
    title (sprintf('LFP Recording, SLE #%d', i));
    ylabel ('mV');
    xlabel ('Time (sec)');
    
%     yyaxis right
%     
%     plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'color', 'k')
%     ylabel ('spike rate/second (Hz)');
%      set(gca,'fontsize',20)
       
    exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
    exportToPPTX('addpicture',figHandle);      
    close(figHandle)
end
        
% save and close the .PPTX
newFile = exportToPPTX('saveandclose',sprintf(excelFileName)); 

end

'successfully completed. Thank you for choosing to use The Epileptiform Detector'
