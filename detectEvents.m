function [ epileptiform, artifacts ] = detectEvents(LFP, minPeakHeight, minPeakDistance)
%UNTITLED4 Summary of this function goes here
%   Input the filtered LFP you want to detect events from


%% Find Light pulse
[P] = pulse_seq(LED);

%Find the quantiles using function quartilesStat
[mx, Q] = quartilesStat(LFP);

%Default values, if minPeakHeight and minPeakDistance is not specified 
if nargin<2
    minPeakHeight = Q(1)*20;   %artifact amplitude >40x 3rd quartile 
    minPeakDistance = 10000;    %artifact spikes seperated by .6 seconds
end

%% Find prominient, distinct spikes in Derivative of filtered LFP (1st search)
[pks_spike, locs_spike] = findpeaks (LFP, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

%% Finding artifacts (Calls function findArtifact.m)
artifacts = findArtifact(LFP, Q(3)*40, 10000);

%% remove artifact spiking (from array of prominient spikes)
for i=1:size(artifacts,1)
    for j=1:numel(locs_spike)
        if locs_spike(j)>=artifacts(i,1) && locs_spike(j)<=artifacts(i,2) 
            locs_spike(j)=-1;
        end
    end
    
end
    %remove spikes that are artifacts
    locs_spike(locs_spike==-1)=[];
       
%% Finding onset 

%Find distance between spikes in data
interSpikeInterval = diff(locs_spike);

%insert "0" into the start interSpikeInterval; allows index to correspond 
n=1;
interSpikeInterval(n+1:end+1,:) = interSpikeInterval(n:end,:);
interSpikeInterval(n,:) = (0);

%Find spikes following 10 s of silence (assume onset)
locs_onset = find (interSpikeInterval(:,1)>100000);

%insert the first epileptiform event into array (not detected with algo)
n=1;
locs_onset(n+1:end+1,:) = locs_onset(n:end,:);
locs_onset(n) = n;
    
%onset times (s)
onsetTimes = zeros(numel (locs_onset),1);
for i=1:numel(locs_onset)
      onsetTimes(i) = t(locs_spike(locs_onset(i)));      
end
               
%% finding Offset 
%should not be a light-triggered spike or an artifact

offsetTimes = zeros(numel (locs_onset),1);

locs_offset = locs_onset - 1;
locs_offset(1) = [];    % there is no spike preceding the very 1st spike

for i=1:numel(locs_offset);
    offsetTimes(i) = t(locs_spike(locs_offset(i)));      
end

%insert the last spike as the offset for last event
offsetTimes(end) = t(locs_spike(end,1));

%% find epileptiform event duration
duration = offsetTimes-onsetTimes;

%putting it all into an array 
epileptiform = [onsetTimes, offsetTimes, duration];

%% Finding light-triggered spikes (events)

%Preallocate
locs_spike(:,2)= 0;

%Find spikes triggered by light
for i=1:numel(P.range(:,1)) %location of light pusle from pulse_seq.m
    range = P.range(i,1):P.range(i,1)+1000; 
    %use or function to combine 2 digital inputs    
    locs_spike(:,2)=or(locs_spike(:,2),ismember (locs_spike(:,1), range));
end

%Store light-triggered spikes
lightTriggeredEvents = locs_spike(locs_spike(:,2)>0, 1);
end

