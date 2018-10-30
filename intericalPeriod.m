function [interictalPeriod, avgDetrendedBaseline, sigmaDetrendedBaseline, figHandle] = intericalPeriod(timeSeries,events,IIS, artifacts, LED, frequency, troubleshooting)
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

timeSeries = LFP_filtered;
interictalPeriod = timeSeries;

%% Indices of interest
epileptiformEventTimes = events(:,1:2);     %Collect all epileptiform events 
epileptiformEventTimes(:,1) = epileptiformEventTimes(:,1) - 1;    %Move onset 0.5s early to make sure all epileptiform activity is accounted for; Warning! error will occur if the first event occured within 0.5 s of recording
epileptiformEventTimes(:,2) = epileptiformEventTimes(:,2) + 3.0;    %Move offset back 3.0s later to make sure all epileptiform activity is accounted for
indexIIEIIS = find(or(events(:,7) == 2, events(:,7) == 3));     %Locate only the IIE & IIS events
epileptiformEventTimes(indexIIEIIS,2) = epileptiformEventTimes(indexIIEIIS,2) + 3.0;  %Move onset back additional 3.0s for IIEs & IISs, the algorithm can't detect their offset effectively
indexFirstSLE = find(events(:,7) == 1, 1, 'first');     %Locate where the first SLE occurs
epileptiformEventTimes = int64(epileptiformEventTimes(indexFirstSLE:end,1:2));     %Ignore all events prior to the first SLE; int64 to make them whole numbers

%% Prepare Time Series 
%Remove spikes (IISs)
for i = 1:size(spikes,1)
    timeStart = int64((spikes(i,1)-1)*frequency);
    timeEnd = int64((spikes(i,2)+6)*frequency);    %Remove 6 s after spike offset
    interictalPeriod(timeStart:timeEnd) = [-1];
    clear timeStart timeEnd
end

%remove artifacts
for i = 1:size(artifactSpikes,1)
    timeStart = int64(artifactSpikes(i,1)*frequency);
    timeEnd = int64(artifactSpikes(i,2)*frequency);    %Remove 6 s after spike offset
    interictalPeriod (timeStart:timeEnd) = [-1];
end

%remove light pulse
if LED
    [pulse] = pulse_seq(LED);   %determine location of light pulses

    %Find range of time when light pulse has potential to trigger an event,
    for i = 1:numel(pulse.range(:,1))
        lightTriggeredOnsetRange = (pulse.range(i,1):pulse.range(i,1)+(6*frequency)); %6 s after light pulse offset
        lightTriggeredOnsetZone{i} = lightTriggeredOnsetRange;
        clear lightTriggeredRange
    end
    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredOnsetZones = cat(2, lightTriggeredOnsetZone{:});  %2 is vertcat

    %% remove spiking due to light pulse
    interictalPeriod (lightTriggeredOnsetZones) = [-1];
end

figure;
reduce_plot(timeSeries)
hold on
reduce_plot(interictalPeriod)

%% Create Vectors of Interictal Period
interictalPeriodCount = numel(epileptiformEventTimes(:,1))-1;   %Period between epileptiform events
interictal = cell(interictalPeriodCount, 1);
for i = 1:interictalPeriodCount
    interictal{i} = interictalPeriod(epileptiformEventTimes(i,2)*frequency:epileptiformEventTimes(i+1,1)*frequency);
    interictal{i} (interictal{i} == -1) = [];   %remove any spikes, artifacfts or like pulses during the interictal period 
    %Characterize baseline features from absolute value of the filtered data
    interictal{i,2} = mean(interictal{i}); %Average
    interictal{i,3} = std(interictal{i}); %Standard Deviation
    figure
    plot (interictal{i})
    title(sprintf('interictal period #%d. Sigma:%.4f', i, interictal{i,3}))
end


% if ~isempty(troubleshooting)
%     figHandle = figure;
%     subplot(2,1,1)
%     reduce_plot(interictalPeriod)
% 
%     subplot(2,1,2)
%     reduce_plot(timeSeries)
% end


end