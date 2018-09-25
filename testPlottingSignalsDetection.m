% %% Data Processing 
% %Center the LFP data
% LFP_centered = LFP - LFP(1);    %Original LFP signal     
% 
% %detrend the LFP data
% LFP_detrended = detrend(LFP); 
% 
% %detrend the LFP data
% LFP_detrendedConstant = detrend(LFP, 'constant'); 
% 
% %Bandpass butter filter [1 - 100 Hz]
% [b,a] = butter(2, [[1 100]/(frequency/2)], 'bandpass');
% LFP_centeredFiltered = filtfilt (b,a,LFP_centered);             %Filtered signal
% 
% %Absolute value of the filtered data
% AbsLFP_centeredFiltered = abs(LFP_centeredFiltered);            %1st derived signal
% 
% %Derivative of the filtered data (absolute value)
% DiffLFP_centeredFiltered = abs(diff(LFP_centeredFiltered));     %2nd derived signal
% 
% %Power of the filtered data (feature for classification)     
% powerFeature = (LFP_centeredFiltered).^2;                     %3rd derived signal
% 

%Putative SLE onset and offset times  
onsetSLE = int64(eventTimes(i,1)*frequency);   %convert into a position
offsetSLE = int64(eventTimes(i,2)*frequency);  %convert into a position

%Baseline adjacent to SLE
onsetBaselineStart = (onsetSLE-(1*frequency));
onsetBaselineEnd = (onsetSLE+(0.5*frequency));
offsetBaselineStart = (offsetSLE-(0.5*frequency));
offsetBaselineEnd = (offsetSLE+(6*frequency));

%Range of LFP to scan for onset
if onsetBaselineStart > 0        
    onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
else
    onsetContext = int64(1:onsetBaselineEnd);
end

%Range of LFP to scan for offset
if offsetBaselineEnd < numel(powerFeatureLowPassFiltered2)
    offsetContext = int64(offsetBaselineStart:offsetBaselineEnd); 
else
    offsetContext = int64(offsetBaselineStart:numel(powerFeatureLowPassFiltered2));
end

%Locating the onset time
prominence = max(powerFeatureLowPassFiltered25(onsetContext))/3; %SLE onset where spike prominience > 1/3 the maximum amplitude
[onset_pks, onset_locs] = findpeaks(powerFeatureLowPassFiltered25(onsetContext), 'MinPeakProminence', prominence);     
if isempty(onset_locs)
    [peakValue, peakIndex] = max(powerFeatureLowPassFiltered25(onsetContext));
    SLEonset_final(i,1) = t(onsetContext(peakIndex)); %The point with maximum power is the onset   
else
    SLEonset_final(i,1) = t(onsetContext(onset_locs(1))); %The onset time, the first spike (increase in power) is the onset   
end

%Locating the offset time    
meanOffsetBaseline = mean(powerFeatureLowPassFiltered2(offsetContext(1:1.5*frequency))); %Mean base line of first 1.5 sec
entireMeanOffsetBaseline = mean(powerFeatureLowPassFiltered2(offsetContext)); %mean baseline of entire signal
OffsetLocation = powerFeatureLowPassFiltered2(offsetContext) > meanOffsetBaseline/2; 
offset_loc = find(OffsetLocation, 1, 'last'); %Last point is the index for the offset location    
offsetSLE_2 = (offsetContext(offset_loc));  %detecting the new offset location         
SLEoffset_final(i,1) = t(offsetSLE_2);

if LED 
    %make sure last spike is not light triggered
    lastBurstEvent = offsetSLE:offsetSLE_2; %the duration of the last burst
    lightTriggered = intersect(lastBurstEvent, lightTriggeredOffsetZones); %check if last burst event is due to a light trigger
    if lightTriggered %if it is light triggered, find last spike that is not light-triggered (new SLE offset)
        %find index of new SLE offset            
        offsetIndex = find(locs_spike(:,1) == int64(eventTimes(i,2)*frequency)); %find the preceding spike, using data from previous 'detectEvent function' in abs LFP
        precedingOffsetIndex = find(flipud(locs_spike(1:offsetIndex, 2))==1,1);        % Thomas flipped the array and looked for the first point where spike was not light-triggered; consider using the while-loop or for-loop backwards
        newOffsetIndex = offsetIndex-(precedingOffsetIndex-1);  %I forgot why I subtract 1, I think it's because it's look at differences
        offsetSLE = int64(locs_spike(newOffsetIndex));  %approximate position of last spike

        %crawl to find the exact offset based on the new offset index
        offsetBaselineStart = double(offsetSLE-(0.5*frequency));  %SLE "context" (start of post-ictal baseline)

        if numel(LFP_filtered)>(offsetSLE+(6*frequency))
            offsetBaselineEnd = double(offsetSLE+(6*frequency));  %SLE "context" (end of post-ictal baseline)
        else
            offsetBaselineEnd = double(numel(powerFeatureLowPassFiltered2)); %SLE "context" (end of post-ictal baseline), end of recording, cut the SLE short
        end

        %Range of LFP to scan 
        offsetContext = (offsetBaselineStart:offsetBaselineEnd);

        %Locating the new offset time    
        meanOffsetBaseline = mean(powerFeatureLowPassFiltered2(offsetContext(1:1.5*frequency))); %Mean base line of first 1.5 sec
        entireMeanOffsetBaseline = mean(powerFeatureLowPassFiltered2(offsetContext)); %mean baseline of entire signal
        offset_loc = find(OffsetLocation, 1, 'last'); %Last point is the offset     
        if offset_loc
            SLEoffset_final(i,1) = t(offsetContext(offset_loc)); %store the detected new offset time             
        else
            SLEoffset_final(i,1) = -1;  %This is considerd to be an artifact and will be removed (revisit this logic), no idea why I did that 
        end
    end
end

%% plotting the onset and offsets detected, troubleshooting purposes     
%% Test plot, onset
%     figure;
%     set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('SLE onset #%d', i)); %select the name you want
%     set(gcf, 'Position', get(0, 'Screensize'));   
%     
%     subplot (3,1,1)
%     plot(t(onsetContext),LFP_filtered(onsetContext))
%     hold on
%     plot(t(onsetSLE), LFP_filtered(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
%     plot(SLEonset_final(i,1), LFP_filtered(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected onset point 
%     plot(t(onsetContext(onset_locs)), LFP_filtered(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     %Labels
%     title ('LFP Bandpass Filtered (1-100 Hz)');
%     ylabel ('mV');
%     xlabel ('Time (sec)');
%     
%     subplot (3,1,2)
%     plot(t(onsetContext),DiffLFP_filtered(onsetContext))
%     hold on
%     plot(t(onsetSLE), DiffLFP_filtered(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
%     plot(SLEonset_final(i,1), DiffLFP_filtered(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected onset point 
%     plot(t(onsetContext(onset_locs)), DiffLFP_filtered(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     %Labels
%     title ('Derivative of LFP Bandpass Filtered (1-100 Hz)');
%     ylabel ('mV');
%     xlabel ('Time (sec)');
%     
%     subplot (3,1,3)
%     plot(t(onsetContext), powerFeatureLowPassFiltered25(onsetContext))
%     hold on
%     plot(t(onsetSLE), powerFeatureLowPassFiltered25(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
%     plot(SLEonset_final(i,1), powerFeatureLowPassFiltered25(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)    %Detected onset point 
%     plot(t(onsetContext(onset_locs)), powerFeatureLowPassFiltered25(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %All potential detected points
%     %Labels
%     title ('Power, Low Pass Filtered (25 Hz)');
%     ylabel ('mV');
%     xlabel ('Time (sec)');    
%     legend ('Low-pass filter, Power', 'Putative Onset', 'detected onsets', 'Final Onset')
%     
%% Test plot, offset
figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name', sprintf ('SLE offset #%d', i)); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));   

subplot (3,1,1)
plot(t(offsetContext),LFP_filtered(offsetContext))
hold on
plot(t(offsetSLE), LFP_filtered(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial putative (rough) detection
plot(SLEoffset_final(i,1), LFP_filtered(offsetContext(offset_loc)), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected offset point 
plot(t(offsetContext(offset_loc)), LFP_filtered(offsetContext(offset_loc)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection

%Labels
title ('LFP Bandpass Filtered (1-100 Hz)');
ylabel ('mV');
xlabel ('Time (sec)');

subplot (3,1,2)
plot(t(offsetContext),DiffLFP_filtered(offsetContext))
hold on
plot(t(offsetSLE), DiffLFP_filtered(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial putative (rough) detection
plot(SLEoffset_final(i,1), DiffLFP_filtered(offsetContext(offset_loc)), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected offset point 
plot(t(offsetContext(offset_loc)), DiffLFP_filtered(offsetContext(offset_loc)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     plot(t(offsetSLE_2), LFP_filtered(offsetSLE_2), '*', 'color', 'blue', 'MarkerSize', 14)  %test 

%Labels
title ('Derivative of LFP Bandpass Filtered (1-100 Hz)');
ylabel ('mV');
xlabel ('Time (sec)');

subplot (3,1,3)
plot(t(offsetContext), powerFeatureLowPassFiltered2(offsetContext))
hold on
plot(t(offsetSLE), powerFeatureLowPassFiltered2(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
plot(SLEoffset_final(i,1), powerFeatureLowPassFiltered2(offsetContext(offset_loc)), 'o', 'color', 'black', 'MarkerSize', 14)    %Detected offset point 
plot(t(offsetContext(offset_loc)), powerFeatureLowPassFiltered2(offsetContext(offset_loc)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection

xL = get(gca, 'XLim');
plot(xL, [meanOffsetBaseline/2 meanOffsetBaseline/2], '--')
plot(xL, [entireMeanOffsetBaseline entireMeanOffsetBaseline], '--')

title ('Power, Low Pass Filtered (2 Hz)');
ylabel ('mV');
xlabel ('Time (sec)');
legend ('Low-pass filter, Power', 'Putative Offset', 'detected offset', 'Final Offset', 'mean baseline/2', 'entire mean baseline')
