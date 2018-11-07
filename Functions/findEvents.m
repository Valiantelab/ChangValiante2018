function [ epileptiformLocation, artifactsLocation, locs_spike ] = findEvents(LFP, frequency, minPeakHeight, minPeakDistance, minArtifactHeight, minArtifactDistance)
%findEvents detects the onset/offset location of events (such as
%epileptiform events, ictal events, interictal spikes, etc.)
%   This function will analyze the time series (i.e., LFP) for major events
%   and artifacts. It is based on the built-in function, findpeaks to
%   detect spikes of a given characteristic. This function detects all the
%   spikes in the time series and assumes that spikes within 10 s of each
%   other are from the same event. Thus, spikes that are seperated by more
%   than 10 seconds represent the space between individual events. The
%   spikes with gaps >10 s represent the onset and offset of neighbouring
%   events. The default criteria is that spikes must have an amplitude that
%   is at least the average value of the time series + 20 x the 1st
%   quartile value, and at least 0.1 seconds apart from each. This allows
%   the algorithm to run faster by reducing the number of spikes detected,
%   but still maintain enough resolution to find the onsets and offsets to
%   the nearest 100 ms. The spikes are considered to be artifacts if they
%   have an amplitude that is greater than 70x the standard deviation
%   (sigma) of the time series above the average value. Detected artifacts 
%   are then ignored when considering which spikes are the onset or 
%   offset of events. Author: Michael Chang (michael.chang@live.ca).

%% Calculate statistics of Time Serise (i.e., LFP recording)
artifact_width = 0.0105;  %seconds 

% Default values, if minPeakHeight and minPeakDistance is not specified 
if nargin<2
    [average, sigma, Q] = statistics(LFP);   %Quartiles
    frequency = 10000;   %10kHz sampling frequency
    minPeakHeight = Q(1)*20;   %spike amplitude >40x 3rd quartile 
    minPeakDistance = 0.1 * frequency;    %spikes seperated by 0.1 seconds
    minArtifactHeight = mean(LFP) + (70*std(LFP));
    minArtifactDistance = 0.6 * frequency;
end

if nargin<3 
    [average, sigma, Q] = statistics(LFP);   %Quartiles Stats
    minPeakHeight = Q(1)*20;   %spike amplitude threshold = 20x(3rd quartile)
    minPeakDistance = 0.1 * frequency;    %spikes seperated by 0.1 seconds
    minArtifactHeight = mean(LFP) + (70*std(LFP));
    minArtifactDistance = 0.6 * frequency;
end


%% Find prominient, distinct spikes in time series of LFP
[pks_spike, locs_spike, width_spike] = findpeaks (LFP, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

if numel(locs_spike) < 2
    fprintf(2,'\nNo epileptiform spikes were detected; review raw data and consider using a lower multiple of baseline sigma as the threshold.\n')
    % Making a empty array so function can complete it's process
    epileptiformLocation = [[], [], []];
    artifactsLocation = [];
else
    %The widths of the spikes may be important when analyzing IIEs and IISs
    locs_spike(:,2) = (width_spike);

    %% Finding artifacts (Calls function findArtifact.m)
    artifactsLocation = findArtifact(LFP, frequency, minArtifactHeight, minArtifactDistance);

    %% remove artifact spiking (from array of prominient spikes)
    for i=1:size(artifactsLocation,1)
        for j=1:numel(locs_spike)
            if locs_spike(j)>=artifactsLocation(i,1) && locs_spike(j)<=artifactsLocation(i,2) 
                locs_spike(j)=-1;
            end
        end    
    end
    %remove spikes that are artifacts
    locs_spike(locs_spike==-1,:)=[];

    %% Finding onset 
    %Find distance between spikes in data
    interSpikeInterval = diff(locs_spike);

    %insert "0" into the start interSpikeInterval; allows index to correspond 
    n=1;
    interSpikeInterval(n+1:end+1,:) = interSpikeInterval(n:end,:);
    interSpikeInterval(n,:) = (0);

    %Find the index of spikes following 10 s of silence (assume onset)
    locs_onset = find (interSpikeInterval(:,1)>100000);

    %insert the first detected spike into array to represent the first event (not detected with algo) 
    n=1;
    locs_onset(n+1:end+1,:) = locs_onset(n:end,:);
    locs_onset(n) = n;   

    %Onset Location(Position)
    onsetLocation = zeros(numel (locs_onset),2);
    for i=1:numel(locs_onset)
          onsetLocation(i,:) = locs_spike(locs_onset(i),:);      
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

    %% Duration
    %Duration of Epileptiform event 
    duration = offsetLocation-onsetLocation(:,1);

    %% Putting onset and offset locations into an array
    epileptiformLocation = [onsetLocation(:,1), offsetLocation, duration, onsetLocation(:,2)];

end

end



