%Program: Epileptiform Activity Detector 
%Author: Michael Chang (michael.chang@live.ca), Fred Chen and Liam Long; 
%Copyright (c) 2018, Valiante Lab
%Version 1.0

%% Clear All
close all
clear all
clc

%% Load .abf and excel data
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', 'F:\');%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
                                                                                         
%% create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%% Seperate signals from .abf files
LFP = x(:,1); 
LED = x(:,2); 

%% normalize the LFP data
LFP_normalized = LFP - LFP(1);

%% Find Light pulse
[P] = pulse_seq(LED);

%% Data Processing 
%Bandpass butter filter [1 - 100 Hz]
[b,a] = butter(2, [[1 100]/(frequency/2)], 'bandpass');
LFP_normalizedFiltered = filtfilt (b,a,LFP_normalized);

%Absolute value of the filtered data 
AbsLFP_normalizedFiltered = abs(LFP_normalizedFiltered);

%Derivative of the filtered data (absolute value)
DiffLFP_normalizedFiltered = abs(diff(LFP_normalizedFiltered));

% %Find the quantiles using function quartilesStat
% [mx, Q] = quartilesStat(DiffLFP_normalizedFiltered);

% %% Find prominient, distinct spikes in Derivative of filtered LFP (1st search)
% [pks_spike, locs_spike] = findpeaks (DiffLFP_normalizedFiltered, 'MinPeakHeight', 20*Q(1), 'MinPeakDistance', 10000);
% 
% %% Find prominient, distinct spikes in (absolute) filtered LFP (2nd search)
% [pks_spike, locs_spike] = findpeaks (AbsLFP_normalizedFiltered, 'MinPeakHeight', Q(3)*1000, 'MinPeakDistance', 1000); 

%% Detect potential events (epileptiform/artifacts) | Derivative Values
[epileptiformLocation, artifacts] = detectEvents (DiffLFP_normalizedFiltered);

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

%Characterize baseline features
sigma = std(AbsLFP_normalizedFilteredBaseline);

%test figure
%figure; plot (AbsLFP_normalizedFilteredBaseline)

%% Detect events (epileptiform/artifacts) | Absolute Values
%Absolute value of the filtered data 
AbsLFP_normalizedFiltered = abs(LFP_normalizedFiltered);

%Detect events
[epileptiformLocation, artifacts] = detectEvents (AbsLFP_normalizedFiltered, 6*sigma, 10000);

%% Finding event time 
%Onset times (s)
onsetTimes = epileptiformLocation(:,1)/frequency; %frequency is your sampling rate

%Offset Times (s)
offsetTimes = epileptiformLocation(:,2)/frequency; 

%Duration of Epileptiform event 
duration = offsetTimes-onsetTimes;

%putting it all into an array 
epileptiform = [onsetTimes, offsetTimes, duration];

%% Identify light-triggered Events
%Find light-triggered spikes 
triggeredSpikes = findTriggeredEvents(AbsLFP_normalizedFiltered, LED);

%Preallocate
epileptiform(:,4)= 0;

%Find light-triggered events 
for i=1:size(epileptiform,1) 
    %use the "ismember" function 
    epileptiform(i,4)=ismember (epileptiformLocation(i,1), triggeredSpikes);
end

%Store light-triggered events (s)
triggeredEvents = epileptiform(epileptiform(:,4)>0, 1);

%% Classifier
SLE = epileptiform(epileptiform(:,3)>=10,:);
IIS = epileptiform(epileptiform(:,3)<10,:);

%% plot graph of normalized  data 

% Find prominient, distinct spikes in (absolute) filtered LFP (2nd search)
[pks_spike, locs_spike] = findpeaks (AbsLFP_normalizedFiltered, 'MinPeakHeight', 6*sigma, 'MinPeakDistance', 1000); 

figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name','Overview of Data'); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

lightpulse = LED > 1;

subplot (3,1,1)
reduce_plot (t, LFP_normalized, 'k');
hold on
reduce_plot (t, lightpulse - 2);

%plot artifacts in red
for i = 1:numel(artifacts(:,1)) 
    reduce_plot (t(artifacts(i,1):artifacts(i,2)), LFP_normalized(artifacts(i,1):artifacts(i,2)), 'r');
end

%plot onset markers
for i=1:numel(epileptiform(:,1))
reduce_plot ((onsetTimes(i)), (LFP_normalized(epileptiformLocation(i))), 'o');
end

%plot offset markers
for i=1:numel(epileptiform(:,2))
reduce_plot ((offsetTimes(i)), (LFP_normalized(epileptiformLocation(i,2))), 'x');
end

title ('Overview of LFP (10000 points/s)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,2) 
reduce_plot (t, AbsLFP_normalizedFiltered, 'b');
hold on
reduce_plot (t, lightpulse - 1);

%plot spikes (will work, with error message; index with artifact removed)
for i=1:size(locs_spike,1)
plot (t(locs_spike(i,1)), (DiffLFP_normalizedFiltered(locs_spike(i,1))), 'x')
end

title ('Overview of filtered LFP (bandpass: 1 to 100 Hz)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,3) 
reduce_plot (t(1:end-1), DiffLFP_normalizedFiltered, 'g');
hold on

% %plot onset markers
% for i=1:numel(locs_onset)
% plot (t(locs_spike(locs_onset(i))), (pks_spike(locs_onset(i))), 'o')
% end
% 
% %plot spikes (will work, with error message; index with artifact removed)
% for i=1:size(locs_spike,1)
% plot (t(locs_spike(i,1)), (DiffLFP_normalizedFiltered(locs_spike(i,1))), 'x')
% end
%  
% %plot offset markers
% for i=1:numel(locs_onset)
% plot ((offsetTimes(i)), (pks_spike(locs_onset(i))), 'x')
% end

title ('Peaks (o) in Derivative of filtered LFP');
ylabel ('Derivative (mV)');
xlabel ('Time (s)');


