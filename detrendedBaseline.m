LFP_detrendedBaseline = LFP_detrended;

%remove events 
for i = 1:size(events,1)
    timeStart = int64(events(i,1)*frequency);
    timeEnd = int64(events (i,2)*frequency);
    LFP_detrendedBaseline(timeStart:timeEnd) = [-1];
    clear timeStart timeEnd
end

%remove IISs
for i = 1:size(IIS,1)
    timeStart = int64(IIS(i,1)*frequency);
    timeEnd = int64(IIS(i,2)*frequency);
    LFP_detrendedBaseline(timeStart:timeEnd) = [-1];
    clear timeStart timeEnd
end

%remove artifacts
for i = 1:size(artifacts,1)
LFP_detrendedBaseline (artifacts(i,1):artifacts(i,2)) = [-1];
end

%remove light pulse
if LED
    [pulse] = pulse_seq(LED);   %determine location of light pulses     

    %Find range of time when light pulse has potential to trigger an event,
    for i = 1:numel(pulse.range(:,1))
        lightTriggeredOnsetRange = (pulse.range(i,1):pulse.range(i,1)+(1*frequency));
        lightTriggeredOnsetZone{i} = lightTriggeredOnsetRange; 
        clear lightTriggeredRange 
    end
    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredOnsetZones = cat(2, lightTriggeredOnsetZone{:});  %2 is vertcat
        
    %% remove spiking due to light pulse 
    LFP_detrendedBaseline (lightTriggeredOnsetZones) = [-1];
end

%Isolate baseline recording
LFP_detrendedBaseline (LFP_detrendedBaseline == -1) = [];

%Characterize baseline features from absolute value of the filtered data 
avgDetrendedBaseline = mean(LFP_detrendedBaseline(2500000:end)); %Average
sigmaDetrendedBaseline = std(LFP_detrendedBaseline(2500000:end)); %Standard Deviation

figure;
reduce_plot(LFP_detrended)
hold on
reduce_plot(LFP_detrendedBaseline)

figure;
subplot(2,1,1)
reduce_plot(LFP_detrendedBaseline(2500000:end))

subplot(2,1,2)
reduce_plot(LFP_detrended(2500000:end))

