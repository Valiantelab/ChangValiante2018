%Program: Epileptiform Activity Detector
%Author: Michael Chang (michael.chang@live.ca)
%Copyright (c) 2018, Valiante Lab
%Version 8.2: Plot distribution of voltage activity's power (workspace) 

%Description: Plots the distributions of voltage activity's power for each
%slice if userInput(5) > 0. Ultimately, it analyzes all the files in the
%folder provided for their ictal period and interictal periods. This script
%will then combine all the periods together and plot their distribution.
%There is a feaible attempt to normalize all the data prior to combining
%them together by dividing the dataset from 1 slice by the smallest voltage
%activity's power. However, I imagine there will be a better way to
%normalize the data in the future that I hope to implement. Final output is
%the figures of the distribution from the ictal periods and interictal
%periods. Future implementations: determine what x_min is and confirm the
%distribution is in fact following a power-law relationship. I imagine that
%x_min will be the border between noise and physiological signals. I need a
%way to objectively define the noise. At the moment, I'm thinking it's fair
%to say noise does not have a dominiant frequency, or perhaps if we can see
%Guassian white noise - that would be the most objective. Note: This script
%will collect interictal periods that contain the interictal events. 

%% Stage 1: Import .Mat Files (workspace)
%clear all (reset)
close all
clear all
clc

%Manually set File Directory
inputdir = 'C:\Users\Michael\OneDrive - University of Toronto\8) Seizure Detection Program\Workspace\Nov2_2018\light-triggered';

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
    [FileName,PathName] = uigetfile ('*.mat','pick .mat file to load Workspace', inputdir);%Choose file    
    fnm = fullfile(PathName,FileName);
    load(sprintf('%s', fnm))   
    userInput(5) = 7;   %troubleshoot
else
    % Analyze all files in folder, multiple files
    PathName = char(InputGUI(6));
    S = dir(fullfile(PathName,'*.mat'));

    for k = 1:numel(S)
        fnm = fullfile(PathName,S(k).name);
        FileName = S(k).name;
        load(sprintf('%s', fnm))     
        userInput(5) = 7;   %troubleshoot
    end
end

%% Stage 2: Process the File
% Author: Michael Chang
% Run this file after the detection algorithm to analyze the results and do
% additional analysis to the detected events. This creats the time vector,
% LFP time series, LED if there is light, and filters the data using a
% bandpass filter (1-50 Hz) and a low pass filter (@68 Hz)

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
%Band Pass Filter
[b,a] = butter(2, ([1 50]/(frequency/2)), 'bandpass');  %Band pass filter
LFP_filteredBandPass = filtfilt (b,a,LFP);             %Bandpass filtered [1 - 50 Hz] singal; because of the 76 Hz noise above, also SLEs only have frequencies up to 20 Hz

%Low Pass Filter
fc = 68; % Cut off frequency
[b,a] = butter(4,fc/(frequency/2), 'low'); %Butterworth filter of order 4
LFP_filtered = filtfilt(b,a,LFP_filteredBandPass); %filtered signal

%Notch Filter
%@60 and @76 Hz

%% Stage 3: Plot Histogram of Voltage Activity's Power 
%Part A: Collect Indices of SLEs and interictal periods
indexSLE = find(events(:,7) == 1);      %Find the indices for all the confirmed SLEs
SLETimes = int64(events(indexSLE,1:2));     %Collect all SLE Times into one array
epileptiformEventTimes = SLETimes;  %Collect all the epileptiform events SLE Times into one array
epileptiformEventTimes(:,1) = epileptiformEventTimes(:,1) - 1;    %Move onset 0.5s early to make sure all epileptiform activity is accounted for; Warning! error will occur if the first event occured within 0.5 s of recording
epileptiformEventTimes(:,2) = epileptiformEventTimes(:,2) + 3.0;    %Move offset back 3.0s later to make sure all epileptiform activity is accounted for

%Part B: Prepare Time Series 
interictalPeriod = LFP_filtered;

%remove artifact Spikes
for i = 1:size(artifactSpikes,1)
    timeStart = int64(artifactSpikes(i,1)*frequency);
    timeEnd = int64(artifactSpikes(i,2)*frequency);    
    interictalPeriod (timeStart:timeEnd) = [-1];
end

%remove artifact Events (events with artifacts)
indexArtifactEvents = find(events(:,7) == 4);
for i = indexArtifactEvents'
    timeStart = int64(events(i,1)*frequency);
    timeEnd = int64((events(i,2)-3)*frequency);    %Remove 3 s after spike offset
    interictalPeriod (timeStart:timeEnd) = [-1];
end

%remove light-induced effects
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

    %remove spiking due to light pulse
    interictalPeriod (lightTriggeredOnsetZones) = [-1];
end

%Part C (1): Create vectors of SLEs 
ictal = cell(numel(SLETimes(:,1)),1);   %preallocate cell array

%Make vectors of SLEs
for i = 1:numel(SLETimes(:,1))
    ictal{i,1} = LFP_filtered(SLETimes(i,1)*frequency:SLETimes(i,2)*frequency);     
end

%Combine together SLE vectors
ictalCombined = ictal(:,1);     %Take the first column from the cell array
ictalCombined = vertcat(ictalCombined{:});  %Vertically concatenate all the cells in the cell array
ictalData{k} = ictalCombined;   %store for Stage 4

%Attempt to normalize combined data | Need to find a better way in the future
data = (abs(ictalCombined));
normalizedIctalData{k} = data/min(data);    %store for Stage 4
clear data

%Part C (2): Create Vectors of Interictal Periods
interictalPeriodCount = numel(epileptiformEventTimes(:,1))-1;   %Period between epileptiform events (period behind last epileptiform event is not 'interictal', technically)
interictal = cell(interictalPeriodCount, 1);    %Preallocate data

%Create Vectors of Interictal Periods
for i = 1:interictalPeriodCount
    interictal{i,1} = interictalPeriod(epileptiformEventTimes(i,2)*frequency:epileptiformEventTimes(i+1,1)*frequency);    %contains IIEs and IISs
    interictal{i,1} (interictal{i,1} == -1) = [];   %remove any artifact (spike/events) or light pulses     
end

%Combine interictal periods together
interictalCombined = interictal(:,1);
interictalCombined = vertcat(interictalCombined {:}); 
interictalData{k} = interictalCombined; 	%store for stage 4

%Attempt to normalize the data | Need to find a better way in the future
data = (abs(interictalCombined));
normalizedInterictalData{k} = data/min(data);       %store for stage 4
clear data

%Part D: Analysis - Plot histogram of distributions of Voltage activity's Power
if userInput(5)>0
    %Plot PowerPoint Slides
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
    exportToPPTX('addtext', 'Distribution of Voltage activity from 4-AP treated cortical slices from VGAT-ChR2', 'Position',[2 1 8 2],...
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
    text = 'The smallest interictal period is 10 s';
    exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 4 5 1],...
                 'Horiz','left', 'Vert','middle', 'FontSize', 14);
    text = 'Data is bandpass filtered 1-100 Hz and Low Pass Filtered with butterworth 4th order filter @68 Hz';
    exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 5 5 1],...
                 'Horiz','left', 'Vert','middle', 'FontSize', 16);

    %% Plot Figures of SLEs
    %Add New Slide
    exportToPPTX('addslide');
    text = 'Distribution of voltage activity from Seizure-Like Events (SLEs)';
    exportToPPTX('addtext', sprintf('%s', text), 'Position',[2 1 8 2],...
                 'Horiz','center', 'Vert','middle', 'FontSize', 36);

    %Plot each SLE
    for i = 1:numel(SLETimes(:,1))    
        figHandle = figure;
        set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
        set(gcf,'Name', sprintf ('SLE #%d', i)); %select the name you want
        set(gcf, 'Position', get(0, 'Screensize'));

        subplot (2,1,1)     %Plot the event vector
        plot (ictal{i})
        title(sprintf('Ictal Event #%d', i))
        axis tight

        subplot (2,1,2)     %Plot the power of the voltage activity
        data = ictal{i,1}.^2;
        histogram(data); %bins the data for you
        set (gca, 'yscale', 'log')
        set (gca, 'xscale', 'log')
        title(sprintf('Ictal Event #%d. Histogram: Distribution of voltage activitys power | Min Data:%.4f  |  Max Data:%.4f ', i, min(data), max(data)))
        xlabel ('Power (mV^2), binned')
        ylabel ('Frequency of Occurrence')

        %Export figures to .pptx
        exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
        exportToPPTX('addpicture',figHandle);
        close(figHandle)
    end

    %Add New Slide
    exportToPPTX('addslide');
    text = 'Distribution of voltage activity from all Seizure-Like Events (SLEs) combined';
    exportToPPTX('addtext', sprintf('%s', text), 'Position',[2 1 8 2],...
                 'Horiz','center', 'Vert','middle', 'FontSize', 36);

    %Plot SLEs Combined
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'SLE combined together'); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));

    subplot (3,1,1)
    plot (ictalCombined)
    title(sprintf('%d Ictal Events combined', numel(ictal(:,1))))
    axis tight

    subplot (3,1,2)
    data = ictalCombined.^2;
    histogram(data);
    set (gca, 'yscale', 'log')
    set (gca, 'xscale', 'log')
    title(sprintf('Histogram: Distribution of Voltage Activity from all Ictal Events. Order of Magnitude difference:%.0fx  |  Min Data:%.4f  |  Max Data:%.4f ', (max(data)/min(data)), min(data), max(data)))
    xlabel ('Power (mV^2), binned')
    ylabel ('Frequency of Occurrence')

    subplot (3,1,3)
    histogram(normalizedIctalData{k}.^2);
    set (gca, 'yscale', 'log')
    set (gca, 'xscale', 'log')
    title(sprintf('Histogram: Normalized Distribution of Voltage Activity from all Ictal Events. Order of Magnitude difference:%.0fx  |  Min Data:%.4f  |  Max Data:%.4f ', (max(data)/min(data)), min(data), max(data)))
    xlabel ('Power (mV^2), binned')
    ylabel ('Frequency of Occurrence')
    
    clear data

    %Export figures to .pptx
    exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
    exportToPPTX('addpicture',figHandle);
    close(figHandle)

    %% Plot Figures of Interictal Periods
    %Add New Slide
    exportToPPTX('addslide');
    text = 'Distribution of voltage activity from Interictal Period (between ictal events)';
    exportToPPTX('addtext', sprintf('%s', text), 'Position',[2 1 8 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 36);
    %Plot Interictal Periods
    for i = 1:interictalPeriodCount
    %Plot figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('Interictal Period #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));

    subplot (2,1,1)
    plot (interictal{i})
    title(sprintf('Interictal Period #%d ', i))
    axis tight

    subplot (2,1,2)
    data = interictal{i}.^2;
    histogram(data);
    set (gca, 'yscale', 'log')
    set (gca, 'xscale', 'log')
    title(sprintf('Interictal Event #%d. Histogram: Distribution of voltage activitys power |  Min Data:%.4f  |  Max Data:%.4f ', i, min(data), max(data)))
    xlabel ('Power (mV^2), binned')
    ylabel ('Frequency of Occurrence')

    %Export figures to .pptx
    exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
    exportToPPTX('addpicture',figHandle);
    close(figHandle)
    end

    %Add New Slide
    exportToPPTX('addslide');
    text = 'Distribution of voltage activity from all Interictal Periods combined';
    exportToPPTX('addtext', sprintf('%s', text), 'Position',[2 1 8 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 36);
    text = 'artifact spikes and events have been removed';
    exportToPPTX('addtext', sprintf('%s', text), 'Position',[3 3 6 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
    exportToPPTX('addtext', 'By: Michael Chang', 'Position',[4 4 4 2],...
                 'Horiz','center', 'Vert','middle', 'FontSize', 20);

    %Plot all interictal periods combined
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Interictal Periods combined together'); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));

    subplot (3,1,1)
    plot (interictalCombined)
    title(sprintf('%d Interictal Periods combined', numel(interictal(:,1))))
    axis tight        

    subplot (3,1,2)
    data = interictalCombined.^2;
    histogram(data);    %this is the power of the voltage activity
    set (gca, 'yscale', 'log')
    set (gca, 'xscale', 'log')
    title(sprintf('Histogram: Distribution of Voltage Activitys power during all interictal periods combined.  Min Data:%.4f  |  Max Data:%.4f ', min(data), max(data)))
    xlabel ('Power (mV^2), binned')
    ylabel ('Frequency of Occurrence')

    subplot (3,1,3)
    data = normalizedInterictalData{k}.^2;
    histogram(data);
    set (gca, 'yscale', 'log')
    set (gca, 'xscale', 'log')
    title(sprintf('Histogram: Normalized Distribution of Voltage Activity from all Ictal Events. Order of Magnitude difference:%.0fx  |  Min Data:%.4f  |  Max Data:%.4f ', (max(data)/min(data)), min(data), max(data)))
    xlabel ('Power (mV^2), binned')
    ylabel ('Frequency of Occurrence')

    %Export figures to .pptx
    exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
    exportToPPTX('addpicture',figHandle);
    close(figHandle)

    clear data
    
    % save and close the .PPTX
    subtitle = '(DistributionVoltagePower)';
    excelFileName = FileName(1:8);
    exportToPPTX('saveandclose',sprintf('%s%s', excelFileName, subtitle));    
end

%end
%end

%% Stage 4: Plot all the data from all recordings 
%Combine ictal periods together
allIctalData = vertcat(ictalData {:});
allNormalizedIctalData = vertcat(normalizedIctalData {:});

%Combine interictal periods together
allInterictalData = vertcat(interictalData {:});
allNormalizedInterictalData = vertcat(normalizedInterictalData {:});

%Plot all distribution combined
  %Plot PowerPoint Slides
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
    exportToPPTX('addtext', 'Distribution of Voltage activity from 4-AP treated cortical slices from VGAT-ChR2', 'Position',[2 1 8 2],...
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
    text = 'The smallest interictal period is 10 s';
    exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 4 5 1],...
                 'Horiz','left', 'Vert','middle', 'FontSize', 14);
    text = 'Data is bandpass filtered 1-100 Hz and Low Pass Filtered with butterworth 4th order filter @68 Hz';
    exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 5 5 1],...
                 'Horiz','left', 'Vert','middle', 'FontSize', 16);

%Add New Slide
exportToPPTX('addslide');
text = 'Distribution of normalized voltage activitys Power from all Ictal Periods combined';
exportToPPTX('addtext', sprintf('%s', text), 'Position',[2 1 8 2],...
         'Horiz','center', 'Vert','middle', 'FontSize', 36);
text = 'artifact spikes and events with artifacts have been removed';
exportToPPTX('addtext', sprintf('%s', text), 'Position',[3 3 6 2],...
         'Horiz','center', 'Vert','middle', 'FontSize', 20);
exportToPPTX('addtext', sprintf('Combined from a total of %d slice recordings', numel(S)), 'Position',[4 4 4 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);

%Plot SLEs Combined
figHandle = figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name', 'ictal periods combined together'); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

% subplot (2,1,1)
% plot (ictalCombined)
% title(sprintf('%d Ictal Events combined', numel(ictal(:,1))))
% axis tight

subplot (2,1,1)
data = allIctalData;
histogram(data.^2);
set (gca, 'yscale', 'log')
set (gca, 'xscale', 'log')
title(sprintf('Histogram: Distribution of raw Voltage Activitys Power from all Ictal Events combined. # of slice recordings :%d  |  Min Data:%.4f  |  Max Data:%.4f ', numel(S), min(data), max(data)))
xlabel ('Power (mV^2), binned')
ylabel ('Frequency of Occurrence')

subplot (2,1,2)
data = allNormalizedIctalData;
histogram(data.^2);
set (gca, 'yscale', 'log')
set (gca, 'xscale', 'log')
title(sprintf('Histogram: Distribution of normalized Voltage Activity Power from all Ictal Events combined. # of slice recordings :%d  |  Min Data:%.4f  |  Max Data:%.4f ', numel(S), min(data), max(data)))
xlabel ('Power (mV^2), binned')
ylabel ('Frequency of Occurrence')


%Export figures to .pptx
exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
exportToPPTX('addpicture',figHandle);
close(figHandle)


%% Interictal Periods Combined
%Add New Slide
exportToPPTX('addslide');
text = 'Distribution of normalized voltage activitys Power from all Interictal Periods combined';
exportToPPTX('addtext', sprintf('%s', text), 'Position',[2 1 8 2],...
         'Horiz','center', 'Vert','middle', 'FontSize', 36);
text = 'artifact spikes and events with artifacts have been removed';
exportToPPTX('addtext', sprintf('%s', text), 'Position',[3 3 6 2],...
         'Horiz','center', 'Vert','middle', 'FontSize', 20);
exportToPPTX('addtext', sprintf('Combined from a total of %d slice recordings', numel(S)), 'Position',[4 4 4 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);

%Plot interictal periods Combined
figHandle = figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name', 'interictal periods combined together'); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

% subplot (2,1,1)
% plot (ictalCombined)
% title(sprintf('%d Ictal Events combined', numel(ictal(:,1))))
% axis tight

subplot (2,1,1)
data = allInterictalData;
histogram(data.^2);
set (gca, 'yscale', 'log')
set (gca, 'xscale', 'log')
title(sprintf('Histogram: Distribution of raw Voltage Activitys Power from all Interictal Periods combined. # of slice recordings :%d  |  Min Data:%.4f  |  Max Data:%.4f ', numel(S), min(data), max(data)))
xlabel ('Power (mV^2), binned')
ylabel ('Frequency of Occurrence')

subplot (2,1,2)
data = allNormalizedInterictalData;
histogram(data.^2);
set (gca, 'yscale', 'log')
set (gca, 'xscale', 'log')
title(sprintf('Histogram: Distribution of normalized Voltage Activitys Power from all Interictal Periods combined. # of slice recordings :%d  |  Min Data:%.4f  |  Max Data:%.4f ', numel(S), min(data), max(data)))
xlabel ('Power (mV^2), binned')
ylabel ('Frequency of Occurrence')

%Export figures to .pptx
exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
exportToPPTX('addpicture',figHandle);
close(figHandle)

% save and close the .PPTX
subtitle = '(CollectiveDistributionVoltagePower)';
excelFileName = 'LightTriggered';
exportToPPTX('saveandclose',sprintf('%s%s', excelFileName, subtitle));    

fprintf(1,'\nThank you for choosing to use the Epileptiform Activity Detector.\n')   
