function [ epileptiformLocation, artifactsLocation, locs_spike ] = detectEvents(LFP, minPeakHeight, minPeakDistance)
%detectEvents is a function to detect the location where events occur
%(i.e., epileptiform events, ictal events, interictal spikes, etc.)
%   This function will take the time series (i.e., LFP) to analyze for
%   major events. It will report the onset and offset of these events based
%   on the default criteria that spikes are 20x 1st quartile in the time
%   series and each spike is seperated by 1 second, unless otherwise
%   specified.

%Find the quantiles using function quartilesStat
[mx, Q] = quartilesStat(LFP);

%Default values, if minPeakHeight and minPeakDistance is not specified 
if nargin<3
    minPeakHeight = Q(1)*20;   %spike amplitude >40x 3rd quartile 
    minPeakDistance = 1000;    %spikes seperated by 0.1 seconds
end

%% Find prominient, distinct spikes in Derivative of filtered LFP (1st search)
[pks_spike, locs_spike] = findpeaks (LFP, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

%% Finding artifacts (Calls function findArtifact.m)
artifactsLocation = findArtifact(LFP);

%% remove artifact spiking (from array of prominient spikes)
for i=1:size(artifactsLocation,1)
    for j=1:numel(locs_spike)
        if locs_spike(j)>=artifactsLocation(i,1) && locs_spike(j)<=artifactsLocation(i,2) 
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

%Onset location
onsetLocation = zeros(numel (locs_onset),1);
for i=1:numel(locs_onset)
      onsetLocation(i) = locs_spike(locs_onset(i));      
end

%% finding Offset 
%should not be a light-triggered spike or an artifact

locs_offset = locs_onset - 1;
locs_offset(1) = [];    % there is no spike preceding the very 1st spike
locs_offset(end+1) = size(locs_spike,1); %insert last spike as last event's offset

%Offset location
offsetLocation = zeros(numel (locs_onset),1);
for i=1:numel(locs_offset);
    offsetLocation(i) = locs_spike(locs_offset(i));      
end

%% Putting onset and offset locations into an array
epileptiformLocation = [onsetLocation, offsetLocation];

end

