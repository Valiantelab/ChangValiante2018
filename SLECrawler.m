function [SLE_final] = SLECrawler(filteredLFP, SLETimes, frequency, LED, onsetDelay, locs_spike, troubleshoot)
%'SLE Crawl' function detects exact onset and offset time of ictal event
%   You upload 1) bandpass filtered LFP data to analyze, 2) the times
%   where all the SLEs (aka ictal events) roughly occur to the nearest 0.5
%   sec, 3) the frequency of the sampling rate. The slecrawl function will
%   then detect the exact onset and offset as Michael Chang would mark the
%   ictal events. Michael determines the valley before the peak (sentinel
%   spike) to be the seizure onset and the point where all spiking activity
%   (power) ends to be the offset. This function can also determine if the
%   SLE is light-triggered and also ensures that the last spike is not  
%   induced by a light pulse. Author: Michael Chang (michael.chang@live.ca)
%   Additional Notes: The default onset delay is 130 ms (or 100 ms after a
%   30 ms light pulse). The threshold for offset detection is
%   'meanOffsetBaseline/2'

   
%% Setting initial variables
%converting inputs into terms Michael used when writing function
LFP_normalizedFiltered = filteredLFP;
SLE = SLETimes;

%Default values, if frequenct is not specified 
if nargin<5
    frequency = 10000;      % 10kHz is default sampling frequency    
    onsetDelay = 0.13;       % seconds after light pusle onset
    troubleshoot = 0;
end

%create time vector
t = (0:(length(filteredLFP)- 1))/frequency;
t = t';

%% Find Light pulse
LED_empty = isempty(LED);
if LED_empty == 0;
    [P] = pulse_seq(LED);

    %determine light pulse onset time
    lightPulseOnset = P.range(:,1)/frequency;

    %Find range of time when light pulse has potential to trigger an event
    %SLE light-triggerd if it starts within 100 ms of light pulse 
    clear lightTriggeredRange lightTriggeredZone
    for i = 1:numel(P.range(:,1))
        lightTriggeredRange = (P.range(i,1):P.range(i,1)+(onsetDelay*frequency));
        lightTriggeredZone{i} = lightTriggeredRange; 
    end

    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredZonesCombined = cat(2, lightTriggeredZone{:});  %2 is vertcat
end
    
%% Processing the data to extract features to determine the onset/offset
%Derivative of the filtered data (absolute value)
DiffLFP_normalizedFiltered = abs(diff(LFP_normalizedFiltered));     %2nd derived signal

%Power of the derivative of the filtered data (absolute values)
powerFeature = (DiffLFP_normalizedFiltered).^2;                     %3rd derived signal

%Lowpass butter filter [2 Hz], to scan for offset
fc = 2; % Cut off frequency
[b,a] = butter(2,fc/(frequency/2)); %Butterworth filter of order 2
powerFeatureLowPassFiltered = filtfilt(b,a,powerFeature); %filtered signal

%Lowpass butter filter [25 Hz], to scan for onset
fc = 25; % Cut off frequency
[b,a] = butter(2,fc/(frequency/2)); %Butterworth filter of order 2
powerFeatureLowPassFiltered25 = filtfilt(b,a,powerFeature); %filtered signal

%% Scanning Low-Pass Filtered Power signal for more accurate onset/offset times
if LED_empty == 1;
    for i = 1:size(SLE,1)
        %Approximate (initially determined) SLE onset and offset times  
        onsetSLE = int64(SLE(i,1)*frequency);   %convert into a position
        offsetSLE = int64(SLE(i,2)*frequency);  %convert into a position

        %SLE "context" (a.k.a. baseline)
        onsetBaselineStart = (onsetSLE-(1*frequency));
        onsetBaselineEnd = (onsetSLE+(0.5*frequency));
        offsetBaselineStart = (offsetSLE-(0.5*frequency));
        offsetBaselineEnd = (offsetSLE+(1*frequency));

        %Range of LFP to scan for onset
        if onsetBaselineStart > 0        
            onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
        else
            onsetContext = int64(1:onsetBaselineEnd);
        end
        
        %Range of LFP to scan for offset
        if offsetBaselineEnd < numel(powerFeatureLowPassFiltered)
            offsetContext = int64(offsetBaselineStart:offsetBaselineEnd); 
        else
            offsetContext = int64(offsetBaselineStart:numel(powerFeatureLowPassFiltered));
        end

        %Locating the onset time
        prominence = max(powerFeatureLowPassFiltered25(onsetContext))/3; %SLE onset where spike prominience > 1/3 the maximum amplitude
        [onset_pks, onset_locs] = findpeaks(powerFeatureLowPassFiltered25(onsetContext), 'MinPeakProminence', prominence);     
        SLEonset_final(i,1) = t(onsetContext(onset_locs(1))); %The onset time, the first spike (increase in power) is the onset   

        %Locating the offset time    
        meanOffsetBaseline = mean(powerFeatureLowPassFiltered(offsetContext)); %SLE ends when signal returned to half the mean power of signal
        OffsetLocation = powerFeatureLowPassFiltered(offsetContext) > meanOffsetBaseline/2; 
        offset_loc = find(OffsetLocation, 1, 'last'); %Last point is the offset     
        offsetSLE_2 = (offsetContext(offset_loc));  %detecting the new offset location         
        SLEoffset_final(i,1) = t(offsetSLE_2);
    end
    
else
    for i = 1:size(SLE,1)
        %Approximate (initially determined) SLE onset and offset times  
        onsetSLE = int64(SLE(i,1)*frequency);   %convert into a position
        offsetSLE = int64(SLE(i,2)*frequency);  %convert into a position

        %SLE "context" (a.k.a. baseline)
        onsetBaselineStart = (onsetSLE-(1*frequency));
        onsetBaselineEnd = (onsetSLE+(0.5*frequency));
        offsetBaselineStart = (offsetSLE-(0.5*frequency));
        offsetBaselineEnd = (offsetSLE+(1*frequency));

        %Range of LFP to scan for onset
        if onsetBaselineStart > 0        
            onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
        else
            onsetContext = int64(1:onsetBaselineEnd);
        end
        
        %Range of LFP to scan for offset
        if offsetBaselineEnd < numel(powerFeatureLowPassFiltered)
            offsetContext = int64(offsetBaselineStart:offsetBaselineEnd); 
        else
            offsetContext = int64(offsetBaselineStart:numel(powerFeatureLowPassFiltered));
        end
        
        %Locating the onset time
        prominence = max(powerFeatureLowPassFiltered25(onsetContext))/3; %SLE onset where spike prominience > 1/3 the maximum amplitude
        [onset_pks, onset_locs] = findpeaks(powerFeatureLowPassFiltered25(onsetContext), 'MinPeakProminence', prominence);     
        if isempty(onset_locs)
            [peakValue, peakIndex] = max(powerFeatureLowPassFiltered25(onsetContext));
            SLEonset_final(i,1) = t(onsetContext(peakIndex)); %The onset time, the first spike (increase in power) is the onset   
        else
            SLEonset_final(i,1) = t(onsetContext(onset_locs(1))); %The onset time, the first spike (increase in power) is the onset   
        end
        
        %Locating the offset time    
        meanOffsetBaseline = mean(powerFeatureLowPassFiltered(offsetContext)); %SLE ends when signal returned to half the mean power of signal
        OffsetLocation = powerFeatureLowPassFiltered(offsetContext) > meanOffsetBaseline/2; 
        offset_loc = find(OffsetLocation, 1, 'last'); %Last point is the offset     
        offsetSLE_2 = (offsetContext(offset_loc));  %detecting the new offset location         
            
        %make sure last spike is not light triggered
        lastBurstEvent = offsetSLE:offsetSLE_2; %the duration of the last burst
        lightTriggered = intersect(lastBurstEvent, lightTriggeredZonesCombined); %check if last burst event is due to a light trigger
        notLightTriggeredEvent = isempty(lightTriggered);
        if notLightTriggeredEvent == 1
            SLEoffset_final(i,1) = t(offsetSLE_2);         
        else   %if it is light triggered find the preceding spike that occured        
            %find index of new SLE offset
            offsetIndex = find(locs_spike == SLE(i,2)*frequency); %find the preceding spike, using data from previous 'detectEvent function' in abs LFP
            
            %find location of new SLE offset
            offsetSLE = int64(locs_spike(offsetIndex-1));  %approximate position of last spike

            %SLE "context" (a.k.a. baseline)
            offsetBaselineStart = (offsetSLE-(0.5*frequency));
            offsetBaselineEnd = (offsetSLE+(1*frequency));

            %Range of LFP to scan 
            offsetContext = int64(offsetBaselineStart:offsetBaselineEnd);

            %Locating the new offset time    
            meanOffsetBaseline = mean(powerFeatureLowPassFiltered(offsetContext)); %SLE ends when signal returned to half the mean power of signal
            OffsetLocation = powerFeatureLowPassFiltered(offsetContext) > meanOffsetBaseline/2; 
            offset_loc = find(OffsetLocation, 1, 'last'); %Last point is the offset     
            SLEoffset_final(i,1) = t(offsetContext(offset_loc)); %store the detect new offset time    
        end
    end

    
    %% plotting the onset and offsets detected, troubleshooting purposes (uncomment to use)
    
    if troubleshoot == 1
    %plot onset
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE onset #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    subplot (2,1,1)
    plot(t(onsetContext),LFP_normalizedFiltered(onsetContext))
    hold on
    plot(t(onsetSLE), LFP_normalizedFiltered(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
    plot(SLEonset_final(i,1), LFP_normalizedFiltered(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected onset point 
    plot(t(onsetContext(onset_locs)), LFP_normalizedFiltered(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('LFP normalized, bandpass filtered');
    ylabel ('mV');
    xlabel ('Time (sec)');
    
    subplot (2,1,2)
    plot(t(onsetContext), powerFeatureLowPassFiltered25(onsetContext))
    hold on
    plot(t(onsetSLE), powerFeatureLowPassFiltered25(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
    plot(SLEonset_final(i,1), powerFeatureLowPassFiltered25(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)    %Detected onset point 
    plot(t(onsetContext(onset_locs)), powerFeatureLowPassFiltered25(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('Power, Low Pass Filtered (2 Hz)');
    ylabel ('mV');
    xlabel ('Time (sec)');
    end
    
end

%Store output 
duration_final = SLEoffset_final - SLEonset_final;
SLE_final = [SLEonset_final, SLEoffset_final, duration_final];  %final list of SLEs, need to filter out artifacts

if LED_empty == 0
    %Preallocate
    SLE_final(:,4)= 0;

    %Classify which SLEs were light triggered
    for i=1:size(SLE_final,1) 
        %use the "ismember" function 
        SLE_final(i,4)=ismember (int64(SLE_final(i,1)*frequency), lightTriggeredZonesCombined);
    end
end



%% Light pulse detect algorithm by Taufik A. Valiante
function [P] = pulse_seq(x)
% Return structure of all the pulses in a sequence of arbitrary height

% Arbitrary threshold
x = x(:);
thresh = max(x)/10;

pulses = x >= thresh;
trans = diff(pulses);
up = find(trans == 1)+1;
down = find(trans == 1)+1;

if isempty(up) || isempty(down)
    P = [];
else
    range = [up,down];
    dur = (down-up);
    isi = [down(1:end-1),up(2:end)];

    P = struct_from_list('range', range, 'dur', dur, 'isi', isi);
end



