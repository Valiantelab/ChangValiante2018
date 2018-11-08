function [timeSeriesBaseline, avgDetrendedBaseline, sigmaDetrendedBaseline, figHandle] = isolateBaseline(timeSeries,events,IIS, artifacts, LED, frequency, troubleshooting)
%baselineIsolation is a function that removes all the events to isolate the
%baseline signal
%   This function removes all the detected events, IISs, Artifacts, to
%   isolate the baseline signal. Specifically, 8 seconds after a
%   epileptiform event, 1 second after a IIS, and 3 second after a light pulse

%Default values if not specified
if nargin < 7
    troubleshooting = [];
    frequency = 1e4;
    LED = [];
end

if nargin < 6
    troubleshooting = [];
    frequency = 1e4;
end

timeSeriesBaseline = timeSeries;

%% Remove Events

%Remove preictal period, preceding 1st ictal event
indexFirstSLE = find(events(:,7) == 1, 1, 'first');
timeStart = 1;
timeEnd = int64(events(indexFirstSLE,1)*frequency);
timeSeriesBaseline(timeStart:timeEnd) = [-1];
clear timeStart timeEnd

%remove events
for i = 1:size(events,1)
    timeStart = int64((events(i,1)-1)*frequency);
    timeEnd = int64((events (i,2)+8)*frequency);
    timeSeriesBaseline(timeStart:timeEnd) = [-1];
    clear timeStart timeEnd
end

%remove IISs
for i = 1:size(IIS,1)
    timeStart = int64((IIS(i,1)-1)*frequency);
    timeEnd = int64((IIS(i,2)+5)*frequency);
    timeSeriesBaseline(timeStart:timeEnd) = [-1];
    clear timeStart timeEnd
end

%remove artifacts
for i = 1:size(artifacts,1)
timeSeriesBaseline (artifacts(i,1):artifacts(i,2)) = [-1];
end

%remove light pulse
if LED
    [pulse] = pulse_seq(LED);   %determine location of light pulses

    %Find range of time when light pulse has potential to trigger an event,
    for i = 1:numel(pulse.range(:,1))
        lightTriggeredOnsetRange = (pulse.range(i,1):pulse.range(i,1)+(5*frequency)); %3 s after light pulse offset
        lightTriggeredOnsetZone{i} = lightTriggeredOnsetRange;
        clear lightTriggeredRange
    end
    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredOnsetZones = cat(2, lightTriggeredOnsetZone{:});  %2 is vertcat

    %% remove spiking due to light pulse
    timeSeriesBaseline (lightTriggeredOnsetZones) = [-1];
end

figure;
reduce_plot(timeSeries)
hold on
reduce_plot(timeSeriesBaseline)

%Isolate baseline recording
timeSeriesBaseline (timeSeriesBaseline == -1) = [];

%Characterize baseline features from absolute value of the filtered data
avgDetrendedBaseline = mean(timeSeriesBaseline); %Average
sigmaDetrendedBaseline = std(timeSeriesBaseline); %Standard Deviation

if ~isempty(troubleshooting)
    figHandle = figure;
    subplot(2,1,1)
    reduce_plot(timeSeriesBaseline)

    subplot(2,1,2)
    reduce_plot(timeSeries)
end


end