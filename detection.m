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

%% Find prominient, distinct spikes in Derivative of filtered LFP
[pks_spike, locs_spike] = findpeaks (DiffLFP_normalizedFiltered, 'MinPeakHeight', 25*Q(1), 'MinPeakDistance', 10000);

%Find distance between spikes in data
interSpikeInterval = diff(locs_spike);

%insert "0" into the start interSpikeInterval, to detect first spike;  
n=1;
interSpikeInterval(n+1:end+1,:) = interSpikeInterval(n:end,:);
interSpikeInterval(n,:) = (0);

%Find spikes following 10 s of silence (assume onset)
[pks_onset, locs_onset] = findpeaks (interSpikeInterval, 'MinPeakHeight', 100000); %Spikes should be at least 10s apart 

%insert first detected spike into onset array
n=1;
locs_onset(n+1:end+1,:) = locs_onset(n:end,:);
locs_onset(n) = n;
    
%% find onset times
onsetTimes = zeros(numel (locs_onset),1);
for i=1:numel(locs_onset)
      onsetTimes(i) = t(locs_spike(locs_onset(i)));      
end

%% find offset times

%find onset time, see if there is another spike for 10 seconds afterwards, 
%that is not light-triggered 
%or an artifact

offsetTimes = zeros(numel (locs_onset),1);

for i=1:numel(locs_onset)
    for j = 1:numel(locs_spike);
        if locs_spike(locs_onset(i)+1) - locs_spike(locs_onset(i)) > 100000
            offsetTimes(i) = locs_onset(i)
        else
            locs_spike(locs_onset(i)+1+j) - locs_spike(locs_onset(i)+j) > 100000
        end
    end
end


for i=1:numel(locs_onset)
    for j = 0:numel(locs_spike);
        if locs_spike(locs_onset(i)+1+j) - locs_spike(locs_onset(i)+j) > 100000
        end
            offsetTimes(i) = locs_onset(i)+j                  
    end
    
end
        


locs_offset = locs_onset - 1;

for i=1:numel(locs_onset);
  
    offsetTimes(i) = t(locs_spike(locs_offset(i)));
       
end

%insert a point in offset array
n=1
insert = t(locs_spike(end))
offsetTimes(end+1) = (insert);

%% find epileptiform event duration
duration = offsetTimes-onsetTimes;

%putting it all into an array 
SLE = [onsetTimes, offsetTimes, duration]

%% Finding artifacts

[pks_artifact, locs_artifact] = findpeaks (DiffLFP_normalizedFiltered, 'MinPeakHeight', Q(3)*30, 'MinPeakDistance', 10000); %artifact should be 30x 3rd quartile 

%preallocate array
artifactStart = zeros(numel(locs_artifact));
artifactEnd = zeros(numel(locs_artifact));
artifactDuration = zeros(numel(locs_artifact));

artifacts = zeros(numel(locs_artifact));



%For-Loop

for i= 1:numel(locs_artifact);

    clear pks_artifact_spikes pks_artifact_spikes
    
    timeSeries=locs_artifact(i)-10000:locs_artifact(i)+10000;

    [pks_artifact_spikes, locs_artifact_spikes] = findpeaks(DiffLFP_normalizedFiltered(timeSeries), 'MinPeakHeight', Q(3)*10); %artifact should be 3x 3rd quartile 

    artifactSpikes=timeSeries(locs_artifact_spikes);

    artifactStart(i) = artifactSpikes(1);
    artifactEnd (i)= artifactSpikes(end);
    artifactDuration(i) = artifactSpikes(end)-artifactSpikes(1);
        
end

    %putting it all into an array 
    artifacts = [artifactStart, artifactEnd, artifactDuration];

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

for i = 1:numel(locs_artifact) %plot artifacts in red
    plot (t(artifacts(i,1):artifacts(i,2)), LFP_normalized(artifacts(i,1):artifacts(i,2)), 'r')
end

%plot onset markers
for i=1:numel(SLE(:,1))
plot (t(locs_spike(locs_onset(i))), (LFP_normalized(locs_onset(i))), 'o')
end

%plot offset markers
for i=1:numel(locs_onset)
plot ((offsetTimes(i)), (LFP_normalized(locs_onset(i))), 'x')
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

%plot onset markers
for i=1:numel(locs_onset)
plot (t(locs_spike(locs_onset(i))), (pks_spike(locs_onset(i))), 'x')
end
 
%plot offset markers
for i=1:numel(locs_onset)
plot ((offsetTimes(i)), (pks_spike(locs_onset(i))), 'o')
end

title ('Peaks (o) in Derivative of filtered LFP');
ylabel ('LFP (mV)');
xlabel ('Time (s)');
