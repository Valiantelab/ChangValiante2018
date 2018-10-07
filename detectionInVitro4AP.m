function [IIS, SLE_final, events] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata)
% inVitro4APDetection is a function designed to search for epileptiform
% events from the in vitro 4-AP seizure model
%   Simply provide the directory to the filename, user inputs, and raw data


%% Standalone function
%Program: Epileptiform Activity Detector 
%Author: Michael Chang (michael.chang@live.ca), Fred Chen and Liam Long; 
%Copyright (c) 2018, Valiante Lab
%Version 6.5

if ~exist('x','var') == 1
    %clear all (reset)
    close all
    clear all
    clc

    %Manually set File DirectorYou seem really sweet and genuine from your profile. 
    inputdir = 'C:\Users\User\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)\1) abf files';

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
    
%% Hard Coded values | Detection settings
%findEvent function
distanceSpike = 0.15;  %distance between spikes (seconds)
distanceArtifact = 0.6; %distance between artifacts (seconds)
minSLEduration = 3.5; %seconds; %change to 5 s if any detection issues 
%SLECrawler function
durationOnsetBaseline = 1.0;     %sec (context to analyze for finding the onset)
durationOffsetBaseline = 1.5;     %sec (context to analyze for finding the offset)
calculateMeanOffsetBaseline = 1.5;    %sec (mean baseline value) | Note: should be smaller than duration

%Label for titles
excelFileName = FileName(1:8);
uniqueTitle = '(detectedEvents)';
finalTitle = '(V6,5)';

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

%% Classify (duration)
%Putative IIS
indexIIS = (epileptiformLocation (:,3)<(minSLEduration*frequency));
epileptiformLocation (indexIIS,7) = 3;   %3 = IIS; 0 = unclassified.

%temp for troublshooting
IIS = epileptiformLocation(indexIIS,1:3)/frequency;

%Putative SLE
indexEvents = (epileptiformLocation (:,3)>=(minSLEduration*frequency));
epileptiformEvents = epileptiformLocation(indexEvents,:);   

%% SLE Crawler: Determine exact onset and offset times | Power Feature
%Scan Low-Pass Filtered Power signal for precise onset/offset times
eventTimes = epileptiformEvents(:,1:2)/frequency;
events = SLECrawler(LFP_filtered, eventTimes, frequency, LED, onsetDelay, offsetDelay, locs_spike_2nd);  

%% Feature Extraction: Duration, Spiking Frequency, Intensity, and Peak-to-Peak Amplitude
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
    
    %Calculate the spiking rate and intensity (per sec) for epileptiform events
    clear spikeRateMinute intensity
    for j = 1:sleDuration
        startWindow = onsetTime+((windowSize*frequency)*(j-1));
        endWindow = onsetTime+((windowSize*frequency)*j);
        %Calculate the spiking rate for epileptiform events
        spikeRate = and(startWindow<=locs_spike_2nd, endWindow >=locs_spike_2nd);
        spikeRateMinute(j,1) = startWindow; %time windows starts
        spikeRateMinute(j,2) = sum(spikeRate(:));   %number of spikes in the window
        %Calculate the intensity per minute for epileptiform events
        if numel(powerFeature) > endWindow
            PowerPerMinute = sum(powerFeature (startWindow:endWindow));        
        else
            PowerPerMinute = sum(powerFeature (startWindow:numel(powerFeature)))
        end                    
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
    
    %Identify epileptiform event phases
    [events(i,13:18), spikeFrequency{i}] = findIctalPhases (spikeFrequency{i});    
    
    %m calculation                   
  
end        

%% Remove mislabelled IIEs and transfer them into IIS array, if present 
indexReclassify = events(:,15) == -1;   %Locate all the IIEs that are less than 3 seconds; please note that this may affect your classification algorithm.

if sum(indexReclassify)>0
    reclassifyEvents = (events(indexReclassify, 1:3));   %Put them into an array
    IIS = [IIS; reclassifyEvents];   %move them into the IIS array
    %Remove IISs from events array
%     events(indexReclassify, :) = [];    %This may affect how your classification algorithm works | may have to revisit
end

%% Classify (extracted features)
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
%indexSLE = events (:,7) == 1;  %Boolean index to indicate which ones are SLEs
indexSLE = find(events (:,7) == 1); %Actual index where SLEs occur, good if you want to perform iterations and maintain position in events
SLE_final = events(indexSLE, [1:8 13]);

%Collect remaining epileptiform event
% indexEpileptiformEvents= events(:,7) == 2; %Boolean index to indicate the remaining epileptiform events
indexEpileptiformEvents= find(events(:,7) == 2); %index to indicate the remaining epileptiform events
epileptiformEvents = events(indexEpileptiformEvents,:);

%% Tertiary Classification, suggested by Taufik to detect Class 5 (questionable) seizures
if epileptiformEvents (:,1) < 1
    fprintf(2,'\nNo questionable epileptiform events were detected in Tertiary Classification.\n')
else if epileptiformEvents(:,1) < 2
         fprintf(2,'\nLimited epileptiform events were detected for Tertiary Classification; epileptiform events will be classified with hard-coded (floor) thresholds.\n')
         epileptiformEvents = classifier_in_vitro (epileptiformEvents, userInput(3));   %Use hardcoded thresholds if there is only 1 event detected       
    else       
        epileptiformEvents = classifier_dynamic (epileptiformEvents, userInput(3), 0.6, 1); %the last input "1" is to indicate it's a IIE classifier
    end
end

%if no Class 5 (questionable) SLEs detected | Repeat classification using hard-coded thresholds, 
if sum(epileptiformEvents (:,7) == 5)<1 || numel(epileptiformEvents (:,7)) < 6    
    fprintf(2,'\nDynamic Classifier did not detect any SLEs. Beginning second attempt with hard-coded thresholds to classify epileptiform events.\n')
    epileptiformEvents = classifier_in_vitro (epileptiformEvents, userInput(3));   %Use hardcoded thresholds if there is only 1 event detected   Also, use in vivo classifier if analying in vivo data        
end

%Recombine outcome into the main events array
events(indexEpileptiformEvents, 7) = epileptiformEvents(:,7);

%% Quaternary Classifier to confirm all SLEs detected have a tonic phase | At Taufik's suggestion
for i = 1:numel(events(:,1))
    %Class 1 SLEs confirmation
    if events(i,7) == 1 && events(i,15) == 0
        events(i,7) = 2;
        fprintf(2,'\nWarning: The Tonic Phase was not detected in SLE #%d from File: %s; so it was reclassified as a IIE.\n', i, FileName)
    
    %IIS confirmation
    else if events(i,7) == 2 && events(i,15) == 0
        events(i,7) = 3;
        fprintf(2,'\nWarning: The Tonic Phase was not detected in Epileptiform Event #%d from File: %s; so it was reclassified as a IIS.\n', i, FileName)
        
        %IIE confirmation
        else if events(i,7) == 5 && events(i,15) == 0
                events(i,7) = 2; 
                 
            else if events(i,7) == 1 && events(i,15) > 0
                    events(i,7) = 1;
                else if events(i,7) == 5 && events(i,15) > 0
                    events(i,7) = 5;                       
                    else
                    events(i,7)= 2;
                    end
                end
            end
        end
    end
end                             

%% Plot all figures for trouble-shooting
if userInput(5) == 1   
    
    %set variables
    timeSeriesLFP = LFP_centered; %Time series to be plotted 

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

    %Plot epileptifrom events detected
    for i = 1:size(events,1)                      
        %Plot Figure
        figHandle = figure;
        set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
        set(gcf,'Name', sprintf ('Epileptiform Event #%d', i)); %select the name you want
        set(gcf, 'Position', get(0, 'Screensize')); 

        subplot (2,1,1)
        figHandle = plotEvent (figHandle, LFP_centered, t, events(i,1:2), locs_spike_2nd, lightpulse);        
        %Labels

        if events(i,7) == 1 
            label = 'SLE';
        else if events(i,7) == 2
                label = 'IIE';
            else if events(i,7) == 3
                    label = 'IIS';
                else if events(i,7) == 4
                        label = 'artifact';
                    else if events(i,15) == -1
                            label = 'IIS';
                            break%revisit
                        else
                            label = 'Questionable SLE';
                        end                    
                    end
                end
            end
        end                   

        title (sprintf('LFP Recording, %s Event #%d | Frequency Change over time', label, i));
        ylabel ('mV');
        xlabel ('Time (sec)');   
        %Plot Frequency Feature
        yyaxis right        

        plot (spikeFrequency{i}(:,1)/frequency, spikeFrequency{i}(:,2), 'o', 'MarkerFaceColor', 'cyan')
  
        if events(i,15) ~= -1    %if it is a IIS, skip this step or it will cause an error
            %set up the index to split frequency feature set
%             activeIndex = spikeFrequency{i}(:,3) == 1;
            inactiveIndex = spikeFrequency{i}(:,3) == 0;
            plot (spikeFrequency{i}(inactiveIndex ,1)/frequency, spikeFrequency{i}(inactiveIndex ,2), 'o', 'MarkerFaceColor', 'magenta')
        end        

        plot ([events(i,13) events(i,13)], ylim)
        plot ([events(i,14) events(i,14)], ylim)

        %Plot again to give markers black border
        plot (spikeFrequency{i}(: ,1)/frequency, spikeFrequency{i}(: ,2), 'o', 'color', 'k')
        
        if events(i,15) == -1    %the extra plot so in the legend we can reveal the classification
            plot (spikeFrequency{i}(: ,1)/frequency, spikeFrequency{i}(: ,2), 'o', 'color', 'k')
        end
        

        %Classification
        if events(i,15) == 0
            classification = 'No tonic phase';
        else if events(i,15) == 1
                classification = 'Tonic-Clonic SLE';
            else if events(i,15) == 2
                    classification = 'Tonic-only SLE';
                else
                    classification = 'IIS';
                end
            end
        end   

        ylabel ('Spike rate/second (Hz)');
        set(gca,'fontsize',14)
        legend ('LFP filtered', 'Epileptiform Event', 'Detected Onset', 'Detected Offset', 'Detected Spikes', 'Applied Stimulus', 'Frequency above threshold', 'Frequency below threshold', 'Tonic Phase onset', 'Clonic Phase onset', sprintf('Classification: %s', classification))
        legend ('Location', 'northeastoutside')

        subplot (2,1,2)
        figHandle = plotEvent (figHandle, LFP_centered, t, events(i,1:2), locs_spike_2nd, lightpulse);                 
        %Labels
        title ('Intensity Change over time');
        ylabel ('mV');
        xlabel ('Time (sec)');   
        %Plot Intensity Feature
        yyaxis right

        %Set up the index to split intensity population
        maxIntensity = double(max(intensityPerMinute{i}(:,2)));
        zeroIndex = intensityPerMinute{i}(:,2) <= (maxIntensity/10); 

        plot (intensityPerMinute{i}(:,1)/frequency, intensityPerMinute{i}(:,2), 'o', 'MarkerFaceColor', 'm')
        plot (intensityPerMinute{i}(zeroIndex ,1)/frequency, intensityPerMinute{i}(zeroIndex ,2), 'o', 'MarkerFaceColor', 'black')
        plot (intensityPerMinute{i}(:,1)/frequency, intensityPerMinute{i}(:,2), 'o', 'color', 'k')      

        ylabel ('intensity/minute');
        set(gca,'fontsize',16)
        set(gca,'fontsize',14)
        legend ('LFP filtered', 'Epileptiform Event', 'Detected Onset', 'Detected Offset', 'Detected Spikes', 'Applied Stimulus', 'High Intensity', 'Low intensity')
        legend ('Location', 'northeastoutside')

        %Export figures to .pptx
        exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
        exportToPPTX('addpicture',figHandle);      
        close(figHandle)

    end

    % save and close the .PPTX
    exportToPPTX('saveandclose',sprintf('%s%s', excelFileName, uniqueTitle)); 
end


%% Collect Class 5 SLEs
indexQuestionableSLE = find(events(:,7)==5);
questionableSLE = events(indexQuestionableSLE,:);

%% IIE Crawler: Determine exact onset and offset times | Power Feature
% indexIIE = events(:,7) == 2;
% IIE_times = events(indexIIE,1:2);
% IIE_final = crawlerIIE(LFP, IIE_times, frequency, LED, onsetDelay, offsetDelay, locs_spike_2nd, 1);  

%% Final Summary Report
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

%% Struct
%File Details
details.FileNameInput = FileName;
details.frequency = frequency;
%User Input and Hard-coded values
details.spikeThreshold = (userInput(1));
details.distanceSpike = distanceSpike;
details.artifactThreshold = userInput(2);
details.distanceArtifact = distanceArtifact;
details.minSLEduration = minSLEduration;
%Detect events (epileptiform/artifacts) | Absolute Values
details.minPeakHeightAbs = minPeakHeight; 
details.minPeakDistanceAbs = minPeakDistance;
details.minArtifactHeightAbs = minArtifactHeight;
details.minArtifactDistanceAbs = minArtifactDistance;
%SLECrawler 
details.durationOnsetBaseline = durationOnsetBaseline;
details.durationOffsetBaseline = durationOffsetBaseline;
details.calculateMeanOffsetBaseline = calculateMeanOffsetBaseline;

%detection results
details.IISsDetected = numel(IIS(:,1));
details.eventsDetected = numel(events(:,1));
details.SLEsDetected = numel (SLE_final(:,1));

%LED details
if userInput(4)>0
    details.stimulusChannel = userInput(4);
    details.onsetDelay = onsetDelay;
    details.offsetDelay = offsetDelay;
end

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

%Sheet 3.5 = Questionable SLE
if isempty(questionableSLE) == 0   
    subtitle4 = {A, B, C, D, E, F, G, H};
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),subtitle4,'QuestionableSLE' ,'A1');
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),questionableSLE(:, 1:8),'QuestionableSLE' ,'A2');
else
    disp ('No Questionable SLEs were detected. Review the raw data to confirm');
end

%Sheet 4 = SLE
if isempty(SLE_final) == 0   
    subtitle4 = {A, B, C, D, E, F, G, H};
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),subtitle4,'SLE' ,'A1');
    xlswrite(sprintf('%s%s',excelFileName, finalTitle ),SLE_final(:, 1:8),'SLE' ,'A2');
else
    disp ('No SLEs were detected. Review the raw data and consider using a lower multiple of baseline sigma as the threshold');
end

%Sheet 0 = Details
writetable(struct2table(details), sprintf('%s%s.xls',excelFileName, finalTitle))    %Plot at end to prevent extra worksheets being produced

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
    
    if LED    
        reduce_plot (t, lightpulse - abs(min(LFP_centered))); %Plot light pulses, if present
    end
    
    plot(xL, [0 0], '--', 'color', [0.5 0.5 0.5])   %Plot dashed line at 0
    
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

    %% Plot figures of classification thresholds
    if userInput(3) == 1    
        figHandles = findall(0, 'Type', 'figure');  %find all open figures exported from classifer
        % Plot figure from Classifer
        for i = numel(figHandles):-1:1      %plot them backwards (in the order they were exported)
            exportToPPTX('addslide'); %Add to new powerpoint slide
            exportToPPTX('addpicture',figHandles(i));      
            close(figHandles(i))
        end    
    end

    %% Plotting out detected SLEs with context | To figure out how off you are    
    for i = 1:size(SLE_final,1) 
       
        %Label figure  
        figHandle = figure;
        set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
        set(gcf,'Name', sprintf ('%s SLE #%d', finalTitle, i)); %select the name you want
        set(gcf, 'Position', get(0, 'Screensize'));  
        
        %Plot figure  
        figHandle = plotEvent (figHandle, LFP_centered, t, SLE_final(i,1:2), locs_spike_2nd, lightpulse); %using Michael's function       
                        
        %Title
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
    exportToPPTX('saveandclose',sprintf('%s(SLEs)', excelFileName)); 
end

fprintf(1,'\nSuccessfully completed. Thank you for choosing to use the In Vitro 4-AP cortical model Epileptiform Detector.\n')



