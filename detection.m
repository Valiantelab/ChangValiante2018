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

%% Find Light pulse
[P] = pulse_seq(LED);

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

%% Finding artifacts

[pks_artifact, locs_artifact] = findpeaks (DiffLFP_normalizedFiltered, 'MinPeakHeight', Q(3)*30, 'MinPeakDistance', 10000); %artifact should be 30x 3rd quartile 

%preallocate array
artifactStart = zeros(size(locs_artifact));
artifactEnd = zeros(size(locs_artifact));
artifactDuration = zeros(size(locs_artifact));


artifacts = zeros(size(locs_artifact,1),3);


%Remove artifact spiking 
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

    %remove artifact spiking from array of prominient spikes
for i=1:size(artifacts,1)
    for j=1:numel(locs_spike)
        if locs_spike(j)>=artifacts(i,1) && locs_spike(j)<=artifacts(i,2) 
            locs_spike(j)=-1;
        end
    end
    
end

    locs_spike(locs_spike==-1)=[];
    
%% Finding onset 

%Find distance between spikes in data
interSpikeInterval = diff(locs_spike);

%insert "0" into the start interSpikeInterval, to detect first spike;  
n=1;
interSpikeInterval(n+1:end+1,:) = interSpikeInterval(n:end,:);
interSpikeInterval(n,:) = (0);

%Find spikes following 10 s of silence (assume onset)
[pks_onset, locs_onset] = findpeaks (interSpikeInterval(:,1), 'MinPeakHeight', 100000); %Spikes should be at least 10s apart 

%insert the first epileptiform event into array (not detected with algo)
n=1;
locs_onset(n+1:end+1,:) = locs_onset(n:end,:);
locs_onset(n) = n;
    
% onset times (s)
onsetTimes = zeros(numel (locs_onset),1);
for i=1:numel(locs_onset)
      onsetTimes(i) = t(locs_spike(locs_onset(i)));      
end

%% Finding light-triggered spikes

%Preallocate
locs_spike(:,2)= 0;

%Find spikes triggered by light
for i=1:numel(P.range(:,1))
    range = P.range(i,1):P.range(i,1)+1000; 
    %use or function to combine 2 digital inputs    
    locs_spike(:,2)=or(locs_spike(:,2),ismember (locs_spike(:,1), range));
end

%Store light-triggered spikes
lightTriggeredEvents = locs_spike(locs_spike(:,2)>0, 1);
                
%% finding Offset 

%find onset time, see if there is another spike for 10 seconds afterwards, 
%that is not light-triggered 
%or an artifact

offsetTimes = zeros(numel (locs_onset),1);

locs_offset = locs_onset - 1;
locs_offset(1) = [];

for i=1:numel(locs_offset);
  
    offsetTimes(i) = t(locs_spike(locs_offset(i)));
       
end

%insert last spike as offset for the last event
insert = t(locs_spike(end,1));
offsetTimes(end) = (insert);


%Fred's stuff (broken)
% 
% for i=1:numel(locs_onset)
%     if (i+1)>numel(locs_onset)
%         interSpikeInterval_sample=interSpikeInterval(locs_onset(i):end);
%         locs_spike_sample=locs_spike(locs_onset(i):end);
%     else
%         interSpikeInterval_sample=interSpikeInterval(locs_onset(i):locs_onset(i+1));
%         locs_spike_sample=locs_spike(locs_onset(i):locs_onset(i+1));
%     end
%     for j=1:numel(interSpikeInterval_sample)
%         currentindex=end-j+1;
%         for k=1:numel(lighttrigger)
%             if(interSpikeInterval_sample (currentindex)-lighttrigger(k))<10000
%                 islastspike=false;
%                 break
%             end
%         end
%         if (islastspike)
%             offsetTimes(i)=locs_spike(currentindex);
%         end
%     end
% end
% 
% 
% for i=1:numel(locs_onset)
%     for j = 0:numel(locs_spike);
%         if locs_spike(locs_onset(i)+1) - locs_spike(locs_onset(i)) > 100000
%             offsetTimes(i) = locs_spike(locs_onset(i)+j)
%         end
%         offsetTimes(i) = locs_onset(i)+j  
%     end
% end
%         



%% find epileptiform event duration
duration = offsetTimes-onsetTimes;

%putting it all into an array 
SLE = [onsetTimes, offsetTimes, duration]

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

%plot artifacts in red
for i = 1:numel(locs_artifact) 
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
plot (t(locs_spike(locs_onset(i))), (pks_spike(locs_onset(i))), 'o')
end

%plot spikes (will work, with error message; index with artifact removed)
for i=1:numel(locs_spike)
plot (t(locs_spike(i)), (DiffLFP_normalizedFiltered(locs_spike(i))), 'x')
end
 
%% test
time_onset=zeros(size(locs_onset))
for i=1:numel(locs_onset)
time=t(locs_spike(locs_onset(i)))
time_onset(i)=time
end

%plot offset markers
for i=1:numel(locs_onset)
plot ((offsetTimes(i)), (pks_spike(locs_onset(i))), 'x')
end


title ('Peaks (o) in Derivative of filtered LFP');
ylabel ('LFP (mV)');
xlabel ('Time (s)');
