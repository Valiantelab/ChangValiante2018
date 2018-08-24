function [ triggeredEvents ] = findTriggeredEvents( LFP, LED, samplingRate, minPeakHeight, minPeakDistance )
%Function: findTriggeredEvents 
%Author: Michael Chang (michael.chang@live.ca)
%Version: 1.0
%Summary: This function finds light triggered events in the
%LFP. It requires the time series for the LFP, the time series of the LED
%(or some other triggered signal), there is an option to specify the height
%of the light triggered event to be considered a spike. 

%Additional details: If the input is not specified, the default values for 
%sampling rate will be 10 kHz,minPeakheight is 20x 1st quartile, and each 
%triggered event must be seperated by 1 second.

%% Delay for triggered onset (seconds)
onsetDelay = .1;

%% Find the quantiles using function quartilesStat
[mx, Q] = quartilesStat(LFP);

%Default values, if minPeakHeight and minPeakDistance is not specified 
if nargin<3
    'Default Values used for finding Triggered events'
    samplingRate = 10000    %smapling rate per second (Frequency)
    minPeakHeight = Q(1)*20   %spike amplitude >40x 3rd quartile 
    minPeakDistance = 1000    %spikes seperated by 1.0 seconds    
end

%% Find Light pulse
[P] = pulse_seq(LED);

%% Find prominient, distinct spikes in Derivative of filtered LFP (1st search)
[pks_spike, locs_spike] = findpeaks (LFP, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

%% Find light-triggered spikes (events)
%Preallocate
locs_spike(:,2)= 0;

%Find light-triggered spikes 
for i=1:numel(P.range(:,1)) %location of light pusle from pulse_seq.m
    range = P.range(i,1):P.range(i,1)+(onsetDelay*samplingRate); 
    %use the "or" function to combine 2 digital inputs    
    locs_spike(:,2)=or(locs_spike(:,2),ismember (locs_spike(:,1), range));
end

%Store light-triggered spikes
triggeredEvents = locs_spike(locs_spike(:,2)>0, 1);
end
