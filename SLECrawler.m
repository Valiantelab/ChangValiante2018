function [SLE_final] = SLECrawler(filteredLFP, SLETimes, frequency, LED, onsetDelay, offsetDelay, locs_spike, troubleshooting)
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
LFP_normalized = filteredLFP;   %for plotting
SLE = SLETimes;

%Default values, if frequenct is not specified 
if nargin<4
    frequency = 10000;      % 10kHz is default sampling frequency    
    onsetDelay = 0.13;      % seconds after light pulse onset to be considered triggered
    offsetDelay = 1.5;      % seconds the event offset follows light pulse to be considered associated
    troubleshooting = 0;    % plot onset and offset detections
end

%create time vector
t = (0:(length(filteredLFP)- 1))/frequency;
t = t';

%% Find Light pulse
if LED;
    [P] = pulse_seq(LED);

    %determine light pulse onset time
    lightPulseOnset = P.range(:,1)/frequency;

    %Find range of time when light pulse has potential to trigger an event
    %SLE light-triggerd if it starts within 100 ms of light pulse 
    clear lightTriggeredRange lightTriggeredZone
    for i = 1:numel(P.range(:,1))
        lightTriggeredOnsetRange = (P.range(i,1):P.range(i,1)+(onsetDelay*frequency));
        lightTriggeredOnsetZone{i} = lightTriggeredOnsetRange; 
        lightTriggeredOffsetRange = (P.range(i,1):P.range(i,1)+(offsetDelay*frequency));
        lightTriggeredOffsetZone{i} = lightTriggeredOffsetRange; 
    end

    %Combine all the ranges where light triggered events occur into one array
    lightTriggeredOnsetZones = cat(2, lightTriggeredOnsetZone{:});  %2 is vertcat
    lightTriggeredOffsetZones = cat(2, lightTriggeredOffsetZone{:});  %2 is vertcat
    
    %% remove spiking due to light pulse (from array of prominient spikes)
    
    locs_spike_replicated = locs_spike(:,1);     %replicate array          
   
    for i = 1:numel(locs_spike_replicated)
        if intersect(locs_spike_replicated(i), lightTriggeredOffsetZones)
            locs_spike_replicated(i)=-1;
        end
    end
        
    %make index of which spikes are due to light pulses
    locs_spike (:,2) = locs_spike_replicated > 0;   %0 means spike triggered by light
    
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
if LED;
   for i = 1:size(SLE,1)
        %Putative SLE onset and offset times  
        onsetSLE = int64(SLE(i,1)*frequency);   %convert into a position
        offsetSLE = int64(SLE(i,2)*frequency);  %convert into a position

        %Baseline adjacent to SLE
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
        lightTriggered = intersect(lastBurstEvent, lightTriggeredOffsetZones); %check if last burst event is due to a light trigger
        if lightTriggered %find location of new SLE offset (that is not light-triggered)
            %find index of new SLE offset            
            offsetIndex = find(locs_spike(:,1) == int64(SLE(i,2)*frequency)); %find the preceding spike, using data from previous 'detectEvent function' in abs LFP
            precedingOffsetIndex = find(flipud(locs_spike(1:offsetIndex, 2))==1,1);        % Thomas flipped the array and looked for the first point where spike was not light-triggered; consider using the while-loop
            newOffsetIndex = offsetIndex-(precedingOffsetIndex-1);
            offsetSLE = int64(locs_spike(newOffsetIndex));  %approximate position of last spike

            %crawl to find the exact offset based on the new offset index
            offsetBaselineStart = double(offsetSLE-(0.5*frequency));  %SLE "context" (preceding baseline)
            
            if numel(filteredLFP)>(offsetSLE+(1*frequency))
                offsetBaselineEnd = double(offsetSLE+(1*frequency));  %SLE "context" (post-ictal baseline)
            else
                offsetBaselineEnd = double(numel(powerFeatureLowPassFiltered)); %SLE "context" (post-ictal baseline), end of recording, cut the SLE short
            end
                       

            %Range of LFP to scan 
            offsetContext = (offsetBaselineStart:offsetBaselineEnd);

            %Locating the new offset time    
            meanOffsetBaseline = mean(powerFeatureLowPassFiltered(offsetContext)); %SLE ends when signal returned to half the mean power of signal
            OffsetLocation = powerFeatureLowPassFiltered(offsetContext) > meanOffsetBaseline/2; 
            offset_loc = find(OffsetLocation, 1, 'last'); %Last point is the offset     
            if offset_loc
                SLEoffset_final(i,1) = t(offsetContext(offset_loc)); %store the detect new offset time             
            else
                SLEoffset_final(i,1) = -1;
            end
            
        else   %if it is not light triggered, use as is 
            SLEoffset_final(i,1) = t(offsetSLE_2);
        end
    end
    

    
    %% plotting the onset and offsets detected, troubleshooting purposes     
    if troubleshooting
    %Test plot, onset
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE onset #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    subplot (2,1,1)
    plot(t(onsetContext),LFP_normalized(onsetContext))
    hold on
    plot(t(onsetSLE), LFP_normalized(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
    plot(SLEonset_final(i,1), LFP_normalized(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected onset point 
    plot(t(onsetContext(onset_locs)), LFP_normalized(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('LFP normalized');
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
    
    %Test plot, offset
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE offset #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    subplot (2,1,1)
    plot(t(offsetContext),LFP_normalized(offsetContext))
    hold on
    plot(t(offsetSLE), LFP_normalized(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial putative (rough) detection
    plot(SLEoffset_final(i,1), LFP_normalized(offsetContext(offset_loc)), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected offset point 
    plot(t(offsetContext(offset_loc)), LFP_normalized(offsetContext(offset_loc)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     plot(t(offsetSLE_2), LFP_normalized(offsetSLE_2), '*', 'color', 'blue', 'MarkerSize', 14)  %test 
    
    %Labels
    title ('LFP normalized');
    ylabel ('mV');
    xlabel ('Time (sec)');
    
    subplot (2,1,2)
    plot(t(offsetContext), powerFeatureLowPassFiltered25(offsetContext))
    hold on
    plot(t(offsetSLE), powerFeatureLowPassFiltered25(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
    plot(SLEoffset_final(i,1), powerFeatureLowPassFiltered25(offsetContext(offset_loc)), 'o', 'color', 'black', 'MarkerSize', 14)    %Detected offset point 
    plot(t(offsetContext(offset_loc)), powerFeatureLowPassFiltered25(offsetContext(offset_loc)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     plot(t(offsetSLE_2), powerFeatureLowPassFiltered25(offsetSLE_2), '*', 'color', 'blue', 'MarkerSize', 14)  %test 
    %Labels
    title ('Power, Low Pass Filtered (2 Hz)');
    ylabel ('mV');
    xlabel ('Time (sec)');
    legend ('Low-pass filter, Power', 'Putative Offset', 'detected offset', 'Final Offset')
    end    
    
else
     for i = 1:size(SLE,1)
        %Putative (initially determined) SLE onset and offset times  
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
end

%Store output 
duration_final = SLEoffset_final - SLEonset_final;
SLE_final = [SLEonset_final, SLEoffset_final, duration_final];  %final list of SLEs, need to filter out artifacts
SLE_final((SLE_final(:,2)==-1),:) = [];     %remove all the rows where SLE is -1

if LED
    %Preallocate
    SLE_final(:,8)= 0;

    %Classify which SLEs were light triggered
    for i=1:size(SLE_final,1) 
        %use the "ismember" function 
        SLE_final(i,8)=ismember (int64(SLE_final(i,1)*frequency), lightTriggeredOnsetZones);
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



