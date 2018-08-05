%Program: Epileptiform Activity Detector 
%Author: Michael Chang (michael.chang@live.ca) and Fred Chen; 
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

%% Data Processing 
%Bandpass butter filter [1 - 100 Hz]
[b,a] = butter(2, [[1 100]/(frequency/2)], 'bandpass');
LFP_normalizedFiltered = filtfilt (b,a,LFP_normalized);

%Derivative of the filtered data (absolute value)
DiffLFP_normalizedFiltered = abs(diff(LFP_normalizedFiltered));

%Find the quantiles using function quartilesStat
[mx, Q] = quartilesStat(DiffLFP_normalizedFiltered);

%Find peaks 
[pks_spike, locs_spike] = findpeaks (DiffLFP_normalizedFiltered, 'MinPeakHeight', 25*Q(1), 'MinPeakDistance', 10000);

%Find start and end of epileptiform events
interSpikeInterval = diff(locs_spike);

%insert a point
n=1;
interSpikeInterval(n+1:end+1,:) = interSpikeInterval(n:end,:);
interSpikeInterval(n,:) = (0);

%Find peaks 
[pks_onset, locs_onset] = findpeaks (interSpikeInterval, 'MinPeakHeight', 100000); %Spikes should be at least 10s apart 

%find onset times
onsetTimes = zeros(numel (locs_onset),1);
for i=1:numel(locs_onset)
  
    onsetTimes(i) = t(locs_spike(locs_onset(i)));
       
 end

%find offset times
offsetTimes = zeros(numel (locs_onset),1);
locs_offset = locs_onset - 1;

for i=1:numel(locs_onset);
  
    offsetTimes(i) = t(locs_spike(locs_offset(i)));
       
end


 %insert a point in onset array
n=1
insert = offsetTimes(1);
onsetTimes(n+1:end+1,:) = onsetTimes(n:end,:);
onsetTimes(n,:) = (insert);

%insert a point in offset array
n=1
insert = t(locs_spike(end))
offsetTimes(end+1) = (insert);

%find epileptiform event duration
duration = offsetTimes-onsetTimes;

%SLE onset and offset times
SLE = [onsetTimes, offsetTimes, duration]

figure
plot (interSpikeInterval);
hold on
plot ((locs_onset), (pks_onset), 'o')



%% plot graph of normalized  data 
figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name','Overview of Data'); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

lightpulse = LED > 1;

subplot (3,1,1)
plot (t, LFP_normalized, 'k')
hold on
plot (t, lightpulse - 2)

%plot onset markers
for i=1:numel(locs_onset)
plot (t(locs_spike(locs_onset(i))), (LFP_normalized(locs_onset(i))), 'x')
end

%plot offset markers
for i=1:numel(locs_onset)
plot ((offsetTimes(i)), (LFP_normalized(locs_onset(i))), 'o')
end

title ('Overview of LFP (10000 points/s)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,2) 
plot (t, LFP_normalizedFiltered, 'b')
hold on
plot (t, lightpulse - 2)
title ('Overview of filtered LFP (bandpass: 1 to 100 Hz)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,3) 
plot (t(1:end-1), DiffLFP_normalizedFiltered, 'g')
hold on
%plot (t(locs_spike), (pks_spike), 'o')
%plot onset markers
for i=1:numel(locs_onset)
plot (t(locs_spike(locs_onset(i))), (pks_spike(locs_onset(i))), 'x')
end
 
title ('Peaks (o) in Derivative of filtered LFP');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

%plot offset markers
for i=1:numel(locs_onset)
plot ((offsetTimes(i)), (pks_spike(locs_onset(i))), 'o')
end
