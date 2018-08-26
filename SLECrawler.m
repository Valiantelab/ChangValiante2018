function [SLE_final] = SLECrawler(filteredLFP, SLETimes,frequency)
%'SLE Crawl' function detects exact onset and offset time of ictal event
%   You upload 1) bandpass filtered LFP data to analyze, 2) the times
%   where all the SLEs (aka ictal events) roughly occur to the nearest 0.5
%   sec, 3) the frequency of the sampling rate. The slecrawl function will
%   then detect the exact onset and offset as Michael Chang would mark the
%   ictal events. Michael determines the valley before the peak (sentinel
%   spike) to be the seizure onset and the point where all spiking activity
%   (power) ends to be the offset. The spike cannot be induced by a light
%   pulse. Author: Michael Chang (michael.chang@live.ca)


%converting inputs into terms Michael used when writing function
LFP_normalizedFiltered = filteredLFP;
SLE = SLETimes;

%Default values, if frequenct is not specified 
if nargin<3
    frequency = 10000;       %10kHz is default sampling frequency    
end

%create time vector
t = (0:(length(filteredLFP)- 1))/frequency;
t = t';

%Processing the data to extract features to determine the onset/offset
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

% Scanning Low-Pass Filtered Power signal for more accurate onset/offset times
for i = 1:size(SLE,1)
    
    %Approximate (initially determined) SLE onset and offset times  
    onsetSLE = int64(SLE(i,1)*frequency);
    offsetSLE = int64(SLE(i,2)*frequency);

    %SLE "context" (a.k.a. baseline)
    onsetBaselineStart = (onsetSLE-(1*frequency));
    onsetBaselineEnd = (onsetSLE+(0.5*frequency));
    offsetBaselineStart = (offsetSLE-(0.5*frequency));
    offsetBaselineEnd = (offsetSLE+(1*frequency));

    %Range of LFP to scan 
    onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
    offsetContext = int64(offsetBaselineStart:offsetBaselineEnd); 

    %Locating the onset time
    prominence = max(powerFeatureLowPassFiltered25(onsetContext))/3; %SLE onset where spike prominience > 1/3 the maximum amplitude
    [onset_pks, onset_locs] = findpeaks(powerFeatureLowPassFiltered25(onsetContext), 'MinPeakProminence', prominence);     
    SLEonset_final(i,1) = t(onsetContext(onset_locs(1))); %The onset time, the first spike (increase in power) is the onset   
   
    %Locating the offset time
    %make sure it's not light triggered
    %[Place Holder]
    
    meanOffsetBaseline = mean(powerFeatureLowPassFiltered(offsetContext)); %SLE ends when signal returned to half the mean power of signal
    OffsetLocation = powerFeatureLowPassFiltered(offsetContext) > meanOffsetBaseline/2; 
    indexOffset = find(OffsetLocation, 1, 'last'); %Last point is the offset     
    SLEoffset_final(i,1) = t(offsetContext(indexOffset)); %the offset time
    

    %% plotting the onset and offsets detected for trouble shooting purposes (uncomment to use)
    
%     %Test plot onset
%     figure;
%     set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('SLE onset #%d', i)); %select the name you want
%     set(gcf, 'Position', get(0, 'Screensize'));   
%     subplot (2,1,1)
%     plot(t(onsetContext),LFP_normalized(onsetContext))
%     hold on
%     plot(t(onsetSLE), LFP_normalized(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
%     plot(SLEonset_final(i,1), LFP_normalized(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected onset point 
%     plot(t(onsetContext(onset_locs)), LFP_normalized(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     %Labels
%     title ('LFP normalized');
%     ylabel ('mV');
%     xlabel ('Time (sec)');
%     
%     subplot (2,1,2)
%     plot(t(onsetContext), powerFeatureLowPassFiltered25(onsetContext))
%     hold on
%     plot(t(onsetSLE), powerFeatureLowPassFiltered25(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
%     plot(SLEonset_final(i,1), powerFeatureLowPassFiltered25(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)    %Detected onset point 
%     plot(t(onsetContext(onset_locs)), powerFeatureLowPassFiltered25(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     %Labels
%     title ('Power, Low Pass Filtered (2 Hz)');
%     ylabel ('mV');
%     xlabel ('Time (sec)');
    
end

duration_final = SLEoffset_final - SLEonset_final;
SLE_final = [SLEonset_final, SLEoffset_final, duration_final];  %final list of SLEs, need to filter out artifacts

end

