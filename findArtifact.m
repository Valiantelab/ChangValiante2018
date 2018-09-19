function [ artifacts, locs_artifact] = findArtifact(LFP, frequency, minPeakHeight, minPeakDistance)

%Program: Artifact Finder
%Authors: Michael Chang (michael.chang@live.ca) and Liam Long
%Description: Searches for artifacts in time series data (bandpass
%     filtered [1 - 100 Hz] LFP. The function is based on the peakfinder
%     function and requires the time series data (LFP), the sampling rate
%     (frequency), the minimum height of the artifact (minPeakheight) and the
%     minimum distance between artifact spikes (minPeakDistance). The default
%     criteria if parameters are not specified is average value of the time
%     series + 70 x the sigma for the height threshold and 0.6 s for the minimum
%     distance between artifact spikes (which tend to be very short in
%     duration). The Output variables will be a matrix containing the details on
%     the location of the artifact onset (1st column); artifact offset (2nd
%     column; and artifact duration (3rd column). The output will be the point
%     in the (not the time) and also the location of the artifcts in terms of
%     position (locs_artifact). To find the time (in seconds), insert the point
%     into the time vectors or simply divide by frequency. A control feature to
%     make sure large biological spikes are not mistaken as a artifact is the
%     maximum width allowed for artifact spikes detected is set to be 10 ms. On
%     average, population spikes tend to be 100 ms in width, while artifacts tend
%     to have a width that is less than 10 ms (on average, ~4 ms). If a
%     putative artifact is detected to be larger than 10.5 ms, it is rejected as a
%     biological spike, and ignored. 


%% Calculate statistics of Time Serise (i.e., LFP recording)
[average, sigma, Q] = statistics(LFP);   %Quartiles, used to detect artifacts on the minor scale

% Default values, if minPeakHeight and minPeakDistance is not specified 
if nargin<2
    frequency = 10000;       %10kHz sampling frequency
    average = mean(LFP);                %Average
    sigma = std(LFP);               %Standard Deviation
    minPeakHeight = average+(sigma*70);   %artifact amplitude >120x 3rd quartile 
    minPeakDistance = 6000;    %artifact spikes seperated by .6 seconds
end

if nargin<3
    average = mean(LFP);                %Average
    sigma = std(LFP);               %Standard Deviation
    minPeakHeight = average+(sigma*70);   %artifact amplitude >120x 3rd quartile 
    minPeakDistance = 6000;    %artifact spikes seperated by .6 seconds
end

%% Finding potential artifacts
[pks_artifact, locs_artifact_potential, width] = findpeaks (LFP, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance); 

%% Finding real Artifacts, width <10 ms
locs_artifact = locs_artifact_potential(width<105); %0.5% lee way provided in threshold

%preallocate array
artifactStart = zeros(size(locs_artifact));
artifactEnd = zeros(size(locs_artifact));
artifactDuration = zeros(size(locs_artifact));
artifacts = zeros(size(locs_artifact,1),3);

%Remove artifact spiking 
for i= 1:numel(locs_artifact);
    clear pks_artifact_spikes pks_artifact_spikes
    
    %Account for artifacts at the very start and end of the recording
     if locs_artifact(i) <= 6000   
         timeSeries = 1:locs_artifact(i)+6000;
     else if (locs_artifact(i)+ 6000) >= numel(LFP)
             timeSeries=locs_artifact(i)-6000:numel(LFP);
         else
             timeSeries=locs_artifact(i)-6000:locs_artifact(i)+6000;
         end
     end        

    [pks_artifact_spikes, locs_artifact_spikes] = findpeaks(LFP(timeSeries), 'MinPeakHeight', Q(3)*5); %artifact should be 3x 3rd quartile 
    
    if locs_artifact_spikes
        artifactSpikes=timeSeries(locs_artifact_spikes);

        artifactStart(i) = artifactSpikes(1);
        artifactEnd (i)= artifactSpikes(end);
        artifactDuration(i) = artifactSpikes(end)-artifactSpikes(1);
    end  
        
end

    %putting it all into an array 
    artifacts = [artifactStart, artifactEnd, artifactDuration];
    
    %test
    
%     figure
%     plot(LFP(timeSeries));
%     hold on
%     plot (locs_artifact_spikes, LFP((locs_artifact_spikes)), '*r')
   

end

