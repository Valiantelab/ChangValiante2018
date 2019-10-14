%% Stage 2: Re-process the File (after manual changes) for further analysis
% Author: Michael Chang
% Run this file after the detection algorithm to analyze the results and do
% additional analysis to the detected events. This creates the time vector,
% LFP time series, LED if there is light, and filters the data using a
% bandpass filter (1-50 Hz) and a low pass filter (@68 Hz)

% % clear all (reset)
% close all
% clear all
% clc

%Add all subfolders in working directory to the path.
addpath(genpath(pwd));  

%Manually set File Directory
inputdir = 'C:\Users\micha\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\4) Acidosis\Mouse 16 - August 8, 2019';

% Load .mat file
    [FileName,PathName] = uigetfile ('*.mat','pick .mat file to load Workspace', inputdir);%Choose file    
    fnm = fullfile(PathName,FileName);
    myVars = {'details', 'samplingInterval', 'x', 'metadata'};
    load(sprintf('%s', fnm), myVars{:})
    S = load(sprintf('%s', fnm));

% f = waitbar(0,'Loading Data: Please wait while data loads...','Name', 'Epileptiform Event Detection in progress');
    
%Load excel file - with human adjusted onset/offset times
    [FileName,PathName] = uigetfile ('*.xls','pick .xls file to load excel sheet', inputdir);%Choose file    
    excel_filename = fullfile(PathName, FileName);
    excelSheet = xlsread(excel_filename, 3);
    %events that were selected after editing the excel sheet by humans
    events = excelSheet(2:end,:);   
    
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
% waitbar(0.02, f, 'Please wait while Processing Data...');

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
    windowSize = 1;  %seconds; you can change window size as desired
    sleDuration = round(numel(eventVector)/frequency);    %rounded to whole number; Note: sometimes the SLE crawler can drop the duration of the event to <1 s
    if sleDuration == 0
        sleDuration = 1;    %need to fix this so you don't analyze event vectors shorter than 1 s
        fprintf(2,'\nWarning! You detected a epileptiform that is shorter than 1 sec, this is an issue with your SLECrawler.m algorithm.')   %testing if I still need this if statement, I think I've fixed teh algorithm so no events <1 s have their features extracted
    end

    %Calculate the intensity (per sec) for epileptiform events; intensity (per sec) is average power
    clear spikeRateMinute intensity    
    for j = 1:sleDuration
        startWindow = onsetTime+((windowSize*frequency)*(j-1));          
        endWindow = onsetTime+((windowSize*frequency)*j);
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
    totalPower(i) = sum(powerFeature(eventVector)); %This is the power of the event
    events (i,5) = totalPower(i) /sleDuration;  %This is the average power/minute

    %Calculate peak-to-peak amplitude of epileptiform event
    eventVectorLFP = LFP_centered(eventVector);
    p2pAmplitude = max(eventVectorLFP) - min (eventVectorLFP);
    events (i,6) = p2pAmplitude;

    %Calculate the event's duration 
    events (i,3) = events (i,2) - events (i,1);
    
    %m calculation
%     test=WP_MultipleRegression(eventVectorLFP', 10000);

    %Calculate dominant frequency
    [epileptiformEvent, interictal] = dominantFrequency(spikes, events, SLE, artifactSpikes, samplingInterval, x, FileName, 2.5, 1.25, 1);

end

%Determine if event is light-triggered, delay to onset
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
    %save final output
    save(sprintf('%s(stage2).mat', FileName(1:end-4)))  %Save Workspace