%function [IIS, SLE_final, events] = inVitro4APDetection(FileName, userInput, x, samplingInterval, metadata)
%inVitro4APDetection is a function designed to search for epileptiform
%events from the in vitro 4-AP seizure model
%   Simply provide the directory to the filename, user inputs, and raw data


%% Standalone function
%Program: Epileptiform Activity Detector 
%Author: Michael Chang (michael.chang@live.ca), Fred Chen and Liam Long; 
%Copyright (c) 2018, Valiante Lab
%Version 6.0

%clear all (reset)
close all
clear all
clc

if ~exist('x','var') == 1
    %Manually set File Director
    inputdir = 'C:\Users\Michael\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)\1) abf files';

    %% GUI to set thresholds
    %Settings, request for user input on threshold
    titleInput = 'Specify Detection Thresholds';
    prompt1 = 'Epileptiform Spike Threshold: average + (3.9 x Sigma)';
    prompt2 = 'Artifact Threshold: average + (70 x Sigma) ';
    prompt3 = 'Figure: Yes (1) or No (0)';
    prompt4 = 'Stimulus channel (enter 0 if none):';
    prompt5 = 'Troubleshooting (plot all epileptiform events): Yes (1) or No 0)';
    prompt6 = 'To analyze multiple files in folder, provide File Directory:';
    prompt = {prompt1, prompt2, prompt3, prompt4, prompt5, prompt6};
    dims = [1 70];
    definput = {'3.9', '70', '0', '2', '0', ''};

    opts = 'on';    %allow end user to resize the GUI window
    InputGUI = (inputdlg(prompt,titleInput,dims,definput, opts));  %GUI to collect End User Inputs
    userInput = str2double(InputGUI(1:5)); %convert inputs into numbers

    %Load .abf file (raw data), analyze single file
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', inputdir);%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
end

    
%% Hard Coded values | distance between spikes
distanceSpike = 0.15;  %distance between spikes (seconds)
distanceArtifact = 0.6; %distance between artifacts (seconds)
minSLEduration = 3.5; %seconds; %change to 5 s if any detection issues 

%Label for titles
excelFileName = FileName(1:8);
uniqueTitle = '(epileptiformEvents)';
finalTitle = '(algoV6)';

%% create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%% Seperate signals from .abf files
LFP = x(:,1);   %original LFP signal
if userInput(4)>0
    LED = x(:,userInput(4));   %light pulse signal, as defined by user's input via GUI
    onsetDelay = 0.13;  %seconds
    offsetDelay = 1.5;  %seconds 
    lightpulse = LED > 1;
else
    LED =[];
    onsetDelay = [];
end

%% Data Processing 
%Center the LFP data
LFP_centered = LFP - LFP(1);                                         

%Bandpass butter filter [1 - 100 Hz]
[b,a] = butter(2, [[1 100]/(frequency/2)], 'bandpass');
LFP_filtered = filtfilt (b,a,LFP);             %Filtered signal

%Absolute value of the filtered data
AbsLFP_Filtered = abs(LFP_filtered);            %1st derived signal

%Derivative of the filtered data (absolute value)
DiffLFP_Filtered = abs(diff(LFP_filtered));     %2nd derived signal

%Power of the filtered data (feature for classification)     
powerFeature = (LFP_filtered).^2;                     %3rd derived signal

%% Detect potential events (epileptiform/artifacts) | Derivative Values
[epileptiformLocation, artifacts, locs_spike_1st] = findEvents (DiffLFP_Filtered, frequency);

%remove potential events
for i = 1:size(epileptiformLocation,1)
AbsLFP_Filtered (epileptiformLocation (i,1):epileptiformLocation (i,2)) = [-1];
end

%remove artifacts
for i = 1:size(artifacts,1)
AbsLFP_Filtered (artifacts(i,1):artifacts(i,2)) = [-1];
end

%Isolate baseline recording
AbsLFP_Filtered (AbsLFP_Filtered == -1) = [];
AbsLFP_centeredFilteredBaseline = AbsLFP_Filtered; %Renamed

%Characterize baseline features from absolute value of the filtered data 
avgBaseline = mean(AbsLFP_centeredFilteredBaseline); %Average
sigmaBaseline = std(AbsLFP_centeredFilteredBaseline); %Standard Deviation

%% Detect events (epileptiform/artifacts) | Absolute Values
%Recreate the Absolute filtered LFP (1st derived signal) vector
AbsLFP_Filtered = abs(LFP_filtered); %the LFP analyzed

%Define thresholds for detection, using inputs from GUI
minPeakHeight = avgBaseline+(userInput(1)*sigmaBaseline);      %threshold for epileptiform spike detection
minPeakDistance = distanceSpike*frequency;                              %minimum distance spikes must be apart
minArtifactHeight = avgBaseline+(userInput(2)*sigmaBaseline);  %threshold for artifact spike detection
minArtifactDistance = distanceArtifact*frequency;                       %minimum distance artifact spikes must be apart

%Detect events
[epileptiformLocation, artifactLocation, locs_spike_2nd] = findEvents (AbsLFP_Filtered, frequency, minPeakHeight, minPeakDistance, minArtifactHeight, minArtifactDistance);

%If no events are detected, terminate script
if isempty(epileptiformLocation)
    fprintf(2,'\nNo epileptiform events were detected. Review the raw data and consider using a different threshold for epileptiform spike detection.\n')
    return
end

%% Stage 0: Initial Classifier (duration)
%Putative IIS
indexIIS = (epileptiformLocation (:,3)<(minSLEduration*frequency));
epileptiformLocation (indexIIS,7) = 3;   %3 = IIS; 0 = unclassified.

%temp for troublshooting
IIS = epileptiformLocation(indexIIS,1:3)/frequency;

%Putative SLE
indexEvents = (epileptiformLocation (:,3)>=(minSLEduration*frequency));
putativeEvents = epileptiformLocation(indexEvents,:);   

%% SLE Crawler: Determine exact onset and offset times | Power Feature
%Scan Low-Pass Filtered Power signal for precise onset/offset times
eventTimes = putativeEvents(:,1:2)/frequency;
events = SLECrawler(LFP_filtered, eventTimes, frequency, LED, onsetDelay, offsetDelay, locs_spike_2nd, 1);  

%% Feature Extraction 
%Optional plot all epileptiform events, for troubleshooting
if userInput(5) == 1   
    
    %set variables
    data1 = LFP_centered; %Time series to be plotted 

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

%Feature Extraction: Duration, Spiking Frequency, Intensity, and Peak-to-peak Amplitude
for i = 1:size(events,1)   
    %make epileptiform event vector
    onsetTime = int64(events(i,1)*frequency);
    offsetTime = int64(events(i,2)*frequency);
    eventVector = int64(onsetTime:offsetTime);  %SLE Vector  
        
    %Split the event vector into (1 min) windows 
    windowSize = 1;  %seconds; you can change window size as desired      
    sleDuration = round(numel(eventVector)/frequency);    %rounded to whole number; Note: sometimes the SLE crawler can drop the duration of the event to <1 s
    if sleDuration == 0
        sleDuration = 1;    %need to fix this so you don't analyze event vectors shorter than 1 s
    end
    
    %Calculate the spiking rate for epileptiform events
    clear spikeRateMinute intensity
    for j = 1:sleDuration
        startWindow = onsetTime+((windowSize*frequency)*(j-1));
        endWindow = onsetTime+((windowSize*frequency)*j);
        %Calculate the spiking rate for epileptiform events
        spikeRate = and(startWindow<=locs_spike_2nd, endWindow >=locs_spike_2nd);
        spikeRateMinute(j,1) = startWindow; %time windows starts
        spikeRateMinute(j,2) = sum(spikeRate(:));   %number of spikes in the window
        %Calculate the intensity per minute for epileptiform events
        PowerPerMinute = sum(powerFeature (startWindow:endWindow));        
        intensity(j,1) = startWindow; %time windows starts
        intensity(j,2) = PowerPerMinute;   %Total power within the (minute) window        
    end
    
    spikeFrequency{i} = spikeRateMinute;    %store the spike frequency of each SLE for plotting later
    intensityPerMinute{i} = intensity;    %store the intensity per minute of each SLE for analysis later
    
    %Calculate average spike rate of epileptiform event
    events (i,4) = mean(spikeRateMinute(:,2));
          
    %Calculate average intensity of epileptiform event
    totalPower(i) = sum(powerFeature(eventVector));
    events (i,5) = totalPower(i) /sleDuration;    
    
    %Calculate peak-to-peak amplitude of epileptiform event
    eventVectorLFP = LFP_centered(eventVector);
    p2pAmplitude = max(eventVectorLFP) - min (eventVectorLFP);
    events (i,6) = p2pAmplitude;
                
    %% Optional: plot vectors for Troubleshooting            
    if userInput(5) == 1    
       
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

        centerLFP = (data1(backgroundVector(1)));
        centerLED = abs(min(data1(backgroundVector)-centerLFP));        
        plot (t(backgroundVector),data1(backgroundVector)-centerLFP)  %background
        hold on
        plot (t(eventVector),data1(eventVector)-centerLFP)     %Epileptiform Event
        plot (t(onsetTime), data1(onsetTime)-centerLFP, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
        plot (t(offsetTime), data1(offsetTime)-centerLFP, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
        indexSpikes = and(onsetTime<locs_spike_2nd, offsetTime>locs_spike_2nd); %Locate spikes between the onset and offset  
        plot (t(locs_spike_2nd(indexSpikes)), (data1(locs_spike_2nd(indexSpikes))-centerLFP), 'x', 'color', 'green') %plot spikes (artifact removed)
        plot (t(backgroundVector),(lightpulse(backgroundVector)/4)-centerLED, 'b') %plot LED 
        
        title (sprintf('LFP Recording, Epileptiform Event #%d | For Troubleshooting', i));
        ylabel ('mV');
        xlabel ('Time (sec)');   

        yyaxis right

%         plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'MarkerFaceColor', 'c')
%         plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'color', 'black')
        plot (intensity(:,1)/frequency, intensity(:,2), 'o', 'MarkerFaceColor', 'c')
        plot (intensity(:,1)/frequency, intensity(:,2), 'o', 'color', 'black')
%         ylabel ('Spike rate/second (Hz)');
        ylabel ('intensity/minute');
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
    
%% Classification  
if events (:,1) < 1
    fprintf(2,'\nNo epileptiform events were detected.\n')
else if events(:,1) < 2
         fprintf(2,'\nLimited Events were detected in the in vitro recording; epileptiform events will be classified with hard-coded thresholds.\n')
         [events, thresholdFrequency, thresholdIntensity, thresholdDuration, indexArtifact, thresholdAmplitudeOutlier] = classifier_in_vitro (events, userInput(3));   %Use hardcoded thresholds if there is only 1 event detected       
    else       
        [events, thresholdFrequency, thresholdIntensity, thresholdDuration, indexArtifact, thresholdAmplitudeOutlier] = classifier_dynamic (events, userInput(3));
    end
end

%if no SLEs detected (because IIEs are absent) | Repeat classification using hard-coded thresholds, 
if sum(events (:,7) == 1)<1 || numel(events (:,7)) < 6    
    fprintf(2,'\nDynamic Classifier did not detect any SLEs. Beginning second attempt with hard-coded thresholds to classify epileptiform events.\n')
    [events, thresholdFrequency, thresholdIntensity, thresholdDuration, indexArtifact, thresholdAmplitudeOutlier] = classifier_in_vitro (events, userInput(3));   %Use hardcoded thresholds if there is only 1 event detected   Also, use in vivo classifier if analying in vivo data        
end
    
%Collect detected SLEs
indexSLE = events (:,7) == 1;  %index to indicate which ones are SLEs
SLE_final = events(indexSLE, :);

%Plot a 3D scatter plot of events
if userInput(3) == 1  
    figEvents = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', '3D Scatter Plot of Epileptiform Events Detected'); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    scatter3(events(:,4), events(:,5), events(:,6), 18, 'black', 'filled')  %All events excluding artifacts
    hold on
    scatter3(events(indexSLE ,4), events(indexSLE ,5), events(indexSLE ,6), 18, 'green', 'filled')   %All SLEs
    scatter3(events(indexArtifact,4), events(indexArtifact,5), events(indexArtifact,6), 18, 'red', 'filled')  %All artifacts    
    %Label
    title ('Classified Epileptiform Events');
    xlabel ('Average Spiking Rate (Hz)');   
    ylabel ('Average Intensity (Power/Duration)');
    zlabel ('Peak-to-Peak Amplitude (mV)');    
    legend ('IIE', 'SLE', 'Artifact')
    legend ('Location', 'southeastoutside')
end
        
% Store light-triggered events (s)
% triggeredEvents = SLE_final(SLE_final(:,4)>0, :);

%% Write to .xls
%set subtitle
A = 'Onset (s)';
B = 'Offset (s)';
C = 'Duration (s)';
D = 'Spike Rate (Hz), average';
E = 'Intensity (mV^2/s), average';
F = 'peak-to-peak Amplitude (mV)';
G = 'Classification';
H = 'Light-triggered';
I = sprintf('%.02f Hz', thresholdFrequency);
J = sprintf('%.02f mV^2/s', thresholdIntensity);
K = sprintf('%.02f s', thresholdDuration);     

%Report outliers detected using peak-to-peak amplitude feature
if sum(indexArtifact)>0
    L = sprintf('%.02f mV', thresholdAmplitudeOutlier);     %artifacts threshold
else
    L = 'no outliers';
end

%Sheet 1 = Details - To be completed at a later date with Liam's help.
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
%     subtitle1 = {details(:,1)};
%     xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),subtitle0,'Details','A1');
%     xlswrite(sprintf('%s%s',excelFileName, uniqueTitle),artifacts/frequency,'Artifacts','A2');

%% Additional Examples of how to display results
% disp(['Mean:                                ',num2str(mx)]);
% disp(['Standard Deviation:                  ',num2str(sigma)]);
% disp(['Median:                              ',num2str(medianx)]);
% disp(['25th Percentile:                     ',num2str(Q(1))]);
% disp(['50th Percentile:                     ',num2str(Q(2))]);
% disp(['75th Percentile:                     ',num2str(Q(3))]);
% disp(['Semi Interquartile Deviation:        ',num2str(SID)]);
% disp(['Number of outliers:                  ',num2str(Noutliers)]);



%Sheet 1 = Events   
if isempty(events) == 0
    subtitle1 = {A, B, C, D, E, F, G, H, I, J, K, L};
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),subtitle1,'Events','A1');
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),events,'Events','A2');
else
    disp ('No events were detected.');
end

%Sheet 2 = Artifacts   
if isempty(artifactLocation) == 0
    subtitle2 = {A, B, C};
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),subtitle2,'Artifacts','A1');
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),artifactLocation/frequency,'Artifacts','A2');
else
    disp ('No artifacts were detected.');
end

%Sheet 3 = IIS
if isempty(IIS) == 0  
    subtitle3 = {A, B, C, G};
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),subtitle3,'IIS' ,'A1');
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),IIS,'IIS' ,'A2');
else
    disp ('No IISs were detected.');
end
    
%Sheet 4 = SLE
if isempty(SLE_final) == 0   
    subtitle4 = {A, B, C, D, E, F, G, H};
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),subtitle4,'SLE' ,'A1');
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),SLE_final(:, 1:8),'SLE' ,'A2');
else
    disp ('No SLEs were detected. Review the raw data and consider using a lower multiple of baseline sigma as the threshold');
end

%% Optional: Plot Figures
if userInput(3) == 1      
%% Creating powerpoint slide
isOpen  = exportToPPTX();
if ~isempty(isOpen)
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

subplot (3,1,1)
reduce_plot (t, LFP_centered, 'k');
hold on
xL = get(gca, 'XLim');  %plot dashed reference line at y = 0
plot(xL, [0 0], '--')

if LED
    reduce_plot (t, lightpulse - abs(min(LFP_centered)));
end

%plot artifacts (red), found in 2nd search
if artifactLocation
    for i = 1:numel(artifactLocation(:,1)) 
        reduce_plot (t(artifactLocation(i,1):artifactLocation(i,2)), LFP_centered(artifactLocation(i,1):artifactLocation(i,2)), 'r');
    end
end

%plot onset/offset markers
if SLE_final(:,1)
    for i=1:numel(SLE_final(:,1))
    reduce_plot ((SLE_final(i,1)), (LFP_centered(int64(SLE_final(i,1)*frequency))), 'o'); %onset markers
    end
    
    for i=1:numel(SLE_final(:,2))
    reduce_plot ((SLE_final(i,2)), (LFP_centered(int64(SLE_final(i,2)*frequency))), 'x'); %offset markers
    end
end

title (sprintf ('Overview of LFP (10000 points/s), %s', FileName));
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,2) 
reduce_plot (t, AbsLFP_Filtered, 'b');
hold on
reduce_plot (t, (lightpulse/2) - 0.75);

%plot spikes (artifact removed)
for i=1:size(locs_spike_2nd,1)
    plot (t(locs_spike_2nd(i,1)), (DiffLFP_Filtered(locs_spike_2nd(i,1))), 'x')
end

title ('Overview of Absolute filtered LFP (bandpass: 1 to 100 Hz)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,3) 
reduce_plot (t(1:end-1), DiffLFP_Filtered, 'g');
hold on

%plot spikes in absoluate derivative of LFP 
for i=1:size(locs_spike_1st,1)
    plot (t(locs_spike_1st(i,1)), (DiffLFP_Filtered(locs_spike_1st(i,1))), 'x')
end

title ('Peaks (o) in Absolute Derivative of filtered LFP');
ylabel ('Derivative (mV)');
xlabel ('Time (s)');

exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
exportToPPTX('addpicture',figHandle);      
close(figHandle)

%% Plot figures of thresholds for classification 
if userInput(3) == 1    
    figHandles = findall(0, 'Type', 'figure');  %find all open figures exported from classifer
    % Plot figure from Classifer
    for i = numel(figHandles):-1:1      %plot them in the order they were exported
        exportToPPTX('addslide'); %Add to new powerpoint slide
        exportToPPTX('addpicture',figHandles(i));      
        close(figHandles(i))
    end    
end

%% Plotting out detected SLEs with context | To figure out how off you are
data1 = LFP_centered; %Time series to be plotted 

for i = 1:size(SLE_final,1) 
    %make SLE vector
    onsetTime = (round(SLE_final(i,1)*10000));
    offsetTime = (round(SLE_final(i,2)*10000));
    sleVector = (onsetTime:offsetTime);  %SLE Vector    
    
    %make background vector
    if (onsetTime >= 50001 && (offsetTime+50000)<numel(data1))
        backgroundVector = int64(onsetTime-50000:offsetTime+50000);   %Background Vector
    elseif (onsetTime < 50001)
        backgroundVector = int64(1:offsetTime+50000);
    elseif ((offsetTime+50000)>numel(data1))
        backgroundVector = int64(onsetTime-50000:numel(data1));
    end
    
    %Plot figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('V4.0 SLE #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));  
    
    centerLFP = (data1(backgroundVector(1)));
    centerLED = abs(min(data1(backgroundVector)-centerLFP));
    plot (t(backgroundVector),data1(backgroundVector)-centerLFP ) %background
    hold on    
    plot (t(sleVector),data1(sleVector)-centerLFP)     %SLE
    plot (t(onsetTime), data1(onsetTime)-centerLFP , 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
    plot (t(offsetTime), data1(offsetTime)-centerLFP , 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
    indexSpikes = and(onsetTime<locs_spike_2nd, offsetTime>locs_spike_2nd); %Locate spikes between the onset and offset  
    plot (t(locs_spike_2nd(indexSpikes)), (data1(locs_spike_2nd(indexSpikes))-centerLFP), 'x', 'color', 'green') %plot spikes (artifact removed)
    if LED
        plot (t(backgroundVector),(lightpulse(backgroundVector))-centerLED , 'b') %plot LED 
    end
       
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
newFile = exportToPPTX('saveandclose',sprintf('%s(SLEs)', excelFileName)); 
end

fprintf(1,'\nSuccessfully completed. Thank you for choosing to use the In Vitro 4-AP cortical model Epileptiform Detector.\n')



