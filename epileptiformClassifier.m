function [ ] = epileptiformClassifier(epileptiformLocation, artifactsLocation, locs_spike, minSLEduration)
%epileptiformClassifier sorts detected epileptiform events into
%seizure-like events (SLEs, aka ictal events), interictal events (IIEs), 
%interictal spikes (IISs), or artifacts.
%   This function is very good because it was written by Michael

%to detect the location where events occur
%(i.e., epileptiform events, ictal events, interictal spikes, etc.)
%   This function will take the time series (i.e., LFP) to analyze for
%   major events. It will report the onset and offset of these events based
%   on the default criteria that spikes are 20x 1st quartile in the time
%   series and each spike is seperated by 1 second, unless otherwise
%   specified.

%% Default values, if frequency and minPeakDistance is not specified 
if nargin<2
    frequency = 10000;   %10kHz sampling frequency
    minSLEduration = 1;  %seconds
end

%convert variables into Michael's terms used when writing "detection.m"
locs_spike_2nd = locs_spike;

%% Feature Extraction (Duration, Spiking Frequency and Intensity)
% Initial Classifier (rough)
putativeSLE = epileptiformLocation(epileptiformLocation(:,3)>=(minSLEduration*frequency),:); 
IIS = epileptiformLocation(epileptiformLocation(:,3)<(minSLEduration*frequency),:); 

for i = 1:size(putativeSLE,1)   
    %make SLE vector
    onsetTime = putativeSLE(i,1);
    offsetTime = putativeSLE(i,2);
    eventVector = (onsetTime:offsetTime);  %SLE Vector  
    
    %Calculate the spiking rate for epileptiform events
    windowSize = 1;  %seconds      
    sleDuration = round(numel(eventVector)/frequency);    %rounded to whole number
    clear spikeRateMinute
    for j = 1:sleDuration
        startWindow = onsetTime+((windowSize*frequency)*(j-1));
        EndWindow = onsetTime+((windowSize*frequency)*j);
        spikeRate = and(startWindow<=locs_spike_2nd, EndWindow >=locs_spike_2nd);
        spikeRateMinute(j,1) = startWindow; %time windows starts
        spikeRateMinute(j,2) = sum(spikeRate(:));   %number of spikes in the window
    end
    
    spikeFrequency{i} = spikeRateMinute;    %store the spike frequency of each SLE for plotting later
    
    %average spike rate of SLE
    putativeSLE (i,4) = mean(spikeRateMinute(:,2));
          
    %average intensity of SLE
    totalPower = sum(powerFeature(eventVector));
    putativeSLE (i,5) = totalPower /sleDuration;     
              
    %make background vector
    if (onsetTime >= 50001 && (offsetTime+50000)<numel(LFP))
        backgroundVector = (onsetTime-50000:offsetTime+50000);   %Background Vector
    elseif (onsetTime < 50001)
        backgroundVector = (1:offsetTime+50000);
    elseif ((offsetTime+50000)>numel(LFP))
        backgroundVector = (onsetTime-50000:numel(LFP));
    end                    

    %% plot vectors
    %use troubleShootingSpikeAnalysis.m
end

%% Classifier - high precision
%Rule #1: average frequency > 1 Hz
index1 = putativeSLE(:,4)>1;

%Rule #2: intensity > (average - sigma)
averageIntensity = mean(putativeSLE(:,5));
sigmaIntensity = std(putativeSLE(:,5));
if averageIntensity > sigmaIntensity 
    thresholdIntensity = averageIntensity - sigmaIntensity; 
else
    thresholdIntensity = averageIntensity; 
end
index2 = putativeSLE(:,5)>thresholdIntensity;

%Rule #3: duration > sigma of durations
sigmaDuration = std(putativeSLE(:,3));
index3 = putativeSLE(:,3)>sigmaDuration; 

%Collect all the SLEs 
SLE = putativeSLE((index1 & index2 & index3), :);   %classified SLEs

%Sort remaining events as interictal events (IIEs)
indexIIS = ~ismember(putativeSLE, SLE);
putativeIIS = putativeSLE(indexIIS(:,1),:); %analyze in future versions


%% finding event time (s)
epileptiformTime = [epileptiformLocation/frequency];

% %% Initial Classifier (rough)
% putativeSLE = epileptiformTime(epileptiformTime(:,3)>=minSLEduration,:); 
% IIS = epileptiformTime(epileptiformTime(:,3)<minSLEduration,:); 

end
