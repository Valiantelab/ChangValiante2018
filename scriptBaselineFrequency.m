%Program: Epileptiform Activity Detector
%Author: Michael Chang (michael.chang@live.ca)
%Copyright (c) 2018, Valiante Lab
%Version 8.1: Locate baseline (Complete)

% Description: Locates segments in the time series that do not have any
% epileptiform activity by looking for sections of actiivty between
% detected events by the algorithm. It will then select the algorithm with
% the lowest sigma to be the best segment of the time series to use as
% baseline. It will analyze the frequency content of the baseline,
% bandpassed between 0-50 Hz; this is because ictal activity in population 
% activity have important frequencies that are below 20 Hz, we don't need
% any of the frequency above that. Furthermore, we have 60 Hz noise and 76
% Hz noise (source unknown) that I don't have to notch filter out anymore.

%% Stage 1: Detect Epileptiform Events
%clear all (reset)
close all
clear all
clc

%Manually set File Directory
inputdir = 'C:\Users\Michael\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)\1) abf files';

%GUI to set thresholds
%Settings, request for user input on threshold
titleInput = 'Specify Detection Thresholds';
prompt1 = 'Epileptiform Spike Threshold: average + (3.9 x Sigma)';
prompt2 = 'Artifact Threshold: average + (70 x Sigma)';
prompt3 = 'Figure: Yes (1) or No (0)';
prompt4 = 'Stimulus channel (enter 0 if none):';
prompt5 = 'Troubleshooting: plot SLEs(1), IIEs(2), IISs(3), Artifacts (4), Review(5), all(6), None(0):';
prompt6 = 'To analyze multiple files in folder, provide File Directory:';
prompt = {prompt1, prompt2, prompt3, prompt4, prompt5, prompt6};
dims = [1 70];
definput = {'3.9', '70', '0', '2', '0', ''};

opts = 'on';    %allow end user to resize the GUI window
InputGUI = (inputdlg(prompt,titleInput,dims,definput, opts));  %GUI to collect End User Inputs
userInput = str2double(InputGUI(1:5)); %convert inputs into numbers

if (InputGUI(6)=="")
    %Load .abf file (raw data), analyze single file
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', inputdir);%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
    [spikes, events, SLE, details, artifactSpikes] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata);
    
    save(sprintf('%s.mat', FileName(1:8)))  %Save Workspace
else
    % Analyze all files in folder, multiple files
    PathName = char(InputGUI(6));
    S = dir(fullfile(PathName,'*.abf'));

    for k = 1:numel(S)
        clear IIS SLE_final events fnm FileName x samplingInterval metadata %clear all the previous data analyzed
        fnm = fullfile(PathName,S(k).name);
        FileName = S(k).name;
        [x,samplingInterval,metadata]=abfload(fnm);
        [spikes, events, SLE, details, artifactSpikes] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata);               
        
        save(sprintf('%s.mat', FileName(1:8)))  %Save Workspace    

%% Stage 2: Process the File
% Author: Michael Chang
% Run this file after the detection algorithm to analyze the results and do
% additional analysis to the detected events. This creats the time vector,
% LFP time series, LED if there is light, and filters the data.

%Create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%Seperate signals from .abf files
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

%Filter Bank
[b,a] = butter(2, ([1 50]/(frequency/2)), 'bandpass');
LFP_filteredBandPass = filtfilt (b,a,LFP);             %Bandpass filtered [1 - 50 Hz] singal; because of the 76 Hz noise above, also SLEs only have frequencies up to 20 Hz

%Low Pass Filter
fc = 68; % Cut off frequency
[b,a] = butter(4,fc/(frequency/2), 'low'); %Butterworth filter of order 4
LFP_filtered = filtfilt(b,a,LFP_filteredBandPass); %filtered signal

%% Stage 3: Find suitable baseline (interictal period with no epileptiform activity)
interictalPeriod = LFP_filtered;    %data analyzed will be the LFP_filtered

%Part A: Indices for interictal periods (between all detected events)
epileptiformEventTimes = events(:,1:2);     %Collect all epileptiform events 
epileptiformEventTimes(:,1) = epileptiformEventTimes(:,1) - 1;    %Move onset 1 s early to make sure all epileptiform activity is accounted for; Warning! error will occur if the first event occured within 0.5 s of recording
epileptiformEventTimes(:,2) = epileptiformEventTimes(:,2) + 3.0;    %Move offset back 3.0s later to make sure all epileptiform activity is accounted for
indexIIEIIS = find(or(events(:,7) == 2, events(:,7) == 3));     %Locate only the IIE & IIS events
epileptiformEventTimes(indexIIEIIS,2) = epileptiformEventTimes(indexIIEIIS,2) + 3.0;  %Move onset back additional 3.0s for IIEs & IISs, the algorithm can't detect their offset effectively
indexFirstSLE = find(events(:,7) == 1, 1, 'first');     %Locate where the first SLE occurs
epileptiformEventTimes = int64(epileptiformEventTimes(indexFirstSLE:end,1:2));     %Ignore all events prior to the first SLE; int64 to make them whole numbers

%Part B: Prepare Time Series 
%Remove spikes (IISs)
for i = 1:size(spikes,1)
    timeStart = int64((spikes(i,1)-1)*frequency);
    timeEnd = int64((spikes(i,2)+6)*frequency);    %Remove 6 s after spike offset
    interictalPeriod(timeStart:timeEnd) = [-1];
    clear timeStart timeEnd
end

%remove artifacts
for i = 1:size(artifactSpikes,1)
    timeStart = int64(artifactSpikes(i,1)*frequency);
    timeEnd = int64(artifactSpikes(i,2)*frequency);    %Remove 6 s after spike offset
    interictalPeriod (timeStart:timeEnd) = [-1];
end

%Note: no need to remove artifact events because they are already accounted
%for in the epileptiformEventTimes

%remove light pulse
if LED
    [pulse] = pulse_seq(LED);   %determine location of light pulses

    %Find range of time when light pulse has potential to trigger an event,
    for i = 1:numel(pulse.range(:,1))
        lightTriggeredOnsetRange = (pulse.range(i,1):pulse.range(i,1)+(6*frequency)); %6 s after light pulse offset
        lightTriggeredOnsetZone{i} = lightTriggeredOnsetRange;
        clear lightTriggeredRange
    end
    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredOnsetZones = cat(2, lightTriggeredOnsetZone{:});  %2 is vertcat

    %% remove spiking due to light pulse
    interictalPeriod (lightTriggeredOnsetZones) = [-1];
end 

%Part C: Create Vectors of Interictal Period
interictalPeriodCount = numel(epileptiformEventTimes(:,1))-1;   %Period between epileptiform events
interictal = cell(interictalPeriodCount, 1);
for i = 1:interictalPeriodCount
    interictal{i} = interictalPeriod(epileptiformEventTimes(i,2)*frequency:epileptiformEventTimes(i+1,1)*frequency);    %Interictal period is between the end of one event and beginning of the next
    interictal{i} (interictal{i} == -1) = [];   %remove any spikes, artifacfts or like pulses during the interictal period 
    if length(interictal{i})<10*frequency
        interictal{i} = -1; %This is a marker to ignore the interictal period <10 s; I only want to analyze periods larger than 10 s
    end
    %Characterize baseline features from absolute value of the filtered data
    interictal{i,2} = mean(interictal{i}); %Average
    interictal{i,3} = std(interictal{i}); %Standard Deviation
%     figure
%     plot (interictal{i})
%     title(sprintf('interictal period #%d. Sigma:%.4f', i, interictal{i,3}))
end

%Locate and delete the interictal period less than 10 secs
indexDelete = find ([interictal{:,2}] == -1); %locate 
interictal(indexDelete,:)=[]; %Delete
clear indexDelete

%Creating powerpoint slide
isOpen  = exportToPPTX();
if ~isempty(isOpen)
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('new','Dimensions',[12 6], ...
    'Title','Epileptiform Event Detector V4.0', ...
    'Author','Michael Chang', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');
%Add New Slide
exportToPPTX('addslide');
exportToPPTX('addtext', 'Frequency Context of Epileptiform Events detected', 'Position',[2 1 8 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 36);
exportToPPTX('addtext', sprintf('File: %s', FileName), 'Position',[3 3 6 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
exportToPPTX('addtext', 'By: Michael Chang', 'Position',[4 4 4 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
%Add New Slide
exportToPPTX('addslide');
exportToPPTX('addtext', 'Legend', 'Position',[0 0 4 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 24);
text = 'Authors: Michael Chang, Liam Long, and Kramay Patel';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 1 6 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = 'Nyquist frequency (Max Frequency/2) typically much higher than physiology frequencies';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 2 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = 'Rayleigh frequency: 1/windowSize (Hz), is the minimum frequency that can be resolved from signal';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 3 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = 'The window size used is 10 s. so the minimum frequncy is 0.1 Hz';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 4 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = 'Accordingly, the smallest event that can be analyzed is 10 s, thus the floor duration for SLE is 10 s';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 5 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 16);

%Part D: Analysis
%Locate the interictal with the lowest sigma, use as baseline
[~, indexMin] = min ([interictal{:,3}]); %locate 

%Plot for your records
i = indexMin;
figHandle = figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf, 'Position', get(0, 'Screensize'));
subplot (2,1,1)
plot (interictal{i})
title(sprintf('Interictal period with lowest Sigma selected to be Baseline | Interictal Period #%d. Sigma: %.4f ', i, interictal{i,3}))
ylabel('Voltage Activity (mV)')
xlabel('data points (Sampling Rate: 10 kHz)')
subplot (2,1,2)
histogram(interictal{i})
title(sprintf('Distribution of voltage activity from Interictal Period #%d', i))
ylabel('Count (Frequency)')
xlabel('Size of Voltage Activity (mV)')

%Export figures to .pptx
exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
exportToPPTX('addpicture',figHandle);
close(figHandle)

%Part 2: Calculate Frequency Content
[nr, ~] = size (interictal);   %Count how many interictal periods there are, "nr"
for i = 1:nr
    %Vector of interictal event with minimual standard deviation 
    eventVector = interictal{i, 1};
    %Energy content of interictal event (serve as baseline to normalize data)
    [s,f,t] = spectrogram (eventVector, 5*frequency, 2.5*frequency, [], frequency, 'yaxis');

    %Dominant Frequency at each time point
    [maxS, idx] = max(abs(s).^2);
    maxFreq = f(idx);

    %decipher
    label = 'Interictal Period (Baseline)';
    if i == indexMin
        classification = sprintf('Minimum Sigma: %.4f (Used as Baseline)',interictal{i,3});
    else
        classification = sprintf('Sigma: %.4f)',interictal{i,3});
    end
    
    %Plot Figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('%s Event #%d', label, i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));

    subplot (3,1,1)
    plot (eventVector)
    title (sprintf('LFP Bandpass Filtered (1-50 Hz), %s Event #%d', label, i))
    xlabel('Data Points')
    ylabel('Voltage (mV)')
    axis tight
    
    subplot (3,2,6)
    histogram (eventVector)
    title (sprintf('Distribution of voltage activity, Bandpass Filtered (1-50 Hz), %s Event #%d', label, i))
    xlabel('Size of Voltage Activity (mV)')
    ylabel('Frequency of Occurance (count)')
    axis tight

    
    subplot (3,1,2)
    contour(t,f,abs(s).^2)
    c = colorbar;
    c.Label.String = 'Power (mV^2)';    %what is the unit really called? 
    ylim([0 100])
    set(gca, 'YScale', 'log')
    title (sprintf('Frequency Content of %s Event #%d. Michaels Algorithm detected: %s', label, i, classification))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')

    subplot (3,2,5)
    plot(t,maxFreq) 
    title (sprintf('Dominant Frequency over duration of %s Event #%d', label, i))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    axis tight
    
     %Export figures to .pptx
     exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
     exportToPPTX('addpicture',figHandle);
     close(figHandle)
end

% save and close the .PPTX
subtitle = '(characterizeBaseline)';
excelFileName = FileName(1:8);
exportToPPTX('saveandclose',sprintf('%s%s', excelFileName, subtitle));

    end
end




fprintf(1,'\nThank you for choosing to use the Valiante Labs Epileptiform Activity Detector.\n')

   