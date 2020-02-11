function [epileptiformEvent, interictal] = dominantFrequency(spikes, events, SLE, artifactSpikes, samplingInterval, x, FileName, windowSize, windowOverlap, figureInput)
%Program: Frequency Content of Epileptiform Events, Work Space 
%Author: Michael Chang (michael.chang@live.ca)
%Copyright (c) 2018, Valiante Lab
%Version 8.2: Analyze Frequency of Epileptiform Events (Complete)

% Description: Function reports the dominant frequency (w/ max PSD) during
% the ictal event. Run function within a script. The data is filted: 
% bandpassed between 0-50 Hz; this is because ictal activity in population 
% activity have important frequencies that are below 20 Hz, we don't need
% any of the frequency above that. Furthermore, we have 60 Hz noise and 76
% Hz noise (source unknown) that I don't have to notch filter out anymore.
% This particular script will open previous workspaces to speed up
% development time.

%% Set parameters for frequency content analysis
if nargin < 6
    windowSize = 2.5;  %seconds; userInput(1)
    windowOverlap = 1.25;   %seconds; userInput(2)
    figureInput = 1; %plot Figure: Yes (1) or No (0)
end
    filter = 'Bandpass: 1-100 Hz'; %Description for subtitles; guiInput{4}
    subtitle = 'frequencyContent';  %set the unique title for .ppt output; guiInput{5}

%% Stage 2: Process the File
% Author: Michael Chang
% Run this file after the detection algorithm to analyze the results and do
% additional analysis to the detected events. This creats the time vector,
% LFP time series, LED if there is light, this stage 2 is unique for
% detecting epileptiform and baseline frequency content

%Create time vector
exactFrequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
frequency = round(exactFrequency);
t = (0:(length(x)- 1))/frequency;
t = t';

%Seperate signals from .abf files
LFP = x(:,1);   %original LFP signal
if size(x,2)>1
    LED = x(:,size(x,2));   %light pulse signal, as defined by user's input via GUI
    onsetDelay = 0.13;  %seconds
    offsetDelay = 1.5;  %seconds
    lightpulse = LED > 1;
else
    LED =[];
    onsetDelay = [];
end

%Filter Bank
% filter = guiInput{4}; %Description for subtitles

%Band Pass Filter
[b,a] = butter(2, ([1 100]/(frequency/2)), 'bandpass');  %Band pass filter
LFP_filteredBandPass = filtfilt (b,a,LFP);             %Bandpass filtered [1 - 100 Hz] singal; because of the 76 Hz noise above, also SLEs only have frequencies up to 20 Hz

%High Pass Filter
fc = 1; % Cut off frequency; a hard stop at 2 Hz
[b,a] = butter(4,fc/(frequency/2), 'high'); %Butterworth filter of order 4
LFP_filteredHighPass = filtfilt(b,a,LFP_filteredBandPass); %filtered signal

%Low Pass Filter
fc = 100; % Cut off frequency; a hard stop at 50 Hz
[b,a] = butter(4,fc/(frequency/2), 'low'); %Bessel filter of order 8
LFP_filteredLowPass = filtfilt(b,a,LFP_filteredHighPass); %filtered signal

%% Stage 3: Find suitable baseline (interictal period with no epileptiform activity), originally used to normalize frequency content 
interictalPeriod = LFP_filteredLowPass;    %data analyzed will be the LFP_filtered

%Part A: Indices for interictal periods (between all detected events)
epileptiformEventTimes = events(:,1:2);     %Collect all epileptiform events 
epileptiformEventTimes(:,1) = epileptiformEventTimes(:,1) - 1;    %Move onset 1 s early to make sure all epileptiform activity is accounted for; Warning! error will occur if the first event occured within 0.5 s of recording
epileptiformEventTimes(:,2) = epileptiformEventTimes(:,2) + 3.0;    %Move offset back 3.0s later to make sure all epileptiform activity is accounted for
indexIIEIIS = find(or(events(:,7) == 2, events(:,7) == 3));     %Locate only the IIE & IIS events
epileptiformEventTimes(indexIIEIIS,2) = epileptiformEventTimes(indexIIEIIS,2) + 3.0;  %Move offset back additional 3.0s for IIEs & IISs, the algorithm can't detect their offset effectively
epileptiformEventTimes = int64(epileptiformEventTimes);     %int64 to make them whole numbers
if ~isempty(SLE) && numel(SLE (:,1)) > 2
    indexFirstSLE = find(events(:,7) == 1, 1, 'first');     %Locate where the first SLE occurs
    epileptiformEventTimes = epileptiformEventTimes(indexFirstSLE:end,1:2);     %Ignore all events prior to the first SLE
else
    disp('No SLEs were detected; Algorithm will calculate interictal period (baseline) from the period between all detected events')
end

%Part B: Prepare Time Series 
%Remove spikes (IISs)
for i = 1:size(spikes,1)
    %Find start time for IIS with padding for content
    if spikes(i,1)> 1
        timeStart = int64((spikes(i,1)-1)*frequency);   %1 second of content before the spike
    else
        timeStart = 1;  %The beginning of the recording, if there's not enough context prior to spike
    end
    %Find end time for IIS with padding for content
    if int64((spikes(i,2)+6)*frequency) < numel(x)
        timeEnd = int64((spikes(i,2)+6)*frequency);    %Remove 6 s after spike offset
    else
        timeEnd = numel(x);    %Take the end of the recording as the end of the IIS
    end
    interictalPeriod(timeStart:timeEnd) = [-1]; %Marked for removal (in part C) 
    clear timeStart timeEnd
end

%remove artifacts
for i = 1:size(artifactSpikes,1)
    timeStart = int64(artifactSpikes(i,1)*frequency);
    timeEnd = int64(artifactSpikes(i,2)*frequency);    %Remove 6 s after spike offset
    interictalPeriod (timeStart:timeEnd) = [-1];    %Marked for removal (in part C)
end

%Note: no need to remove artifact events because they are already accounted for in the epileptiformEventTimes

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

    %remove spiking due to light pulse
    interictalPeriod (lightTriggeredOnsetZones) = [-1]; %Marked for removal 
end 

%Part C (1): Create vectors of detected events 
% indexEvent = find(events(:,3)>10);
% eventTimes = int64(events(indexEvent,:));
epileptiformEvent = cell(numel(events(:,1)),3);   %preallocate cell array

%Pad the epileptiform event vector equal to the windowSize
% windowSize = userInput(1);
% windowOverlap = userInput(2);
minInterictalPeriod = windowSize*2;   %secs

%Make vectors of SLEs 
for i = 1:numel(events(:,1))  
    [~, indicesBackground] = eventIndices(LFP_filteredLowPass, events(i,:), windowSize, frequency);    %Make vectors based on original times detected by algorithm
    epileptiformEvent{i,1} = LFP_filteredLowPass(indicesBackground); 
end

%Part C(2): Create Vectors of Interictal Period
interictalPeriodCount = numel(epileptiformEventTimes(:,1))-1;   %Period between epileptiform events
interictal = cell(interictalPeriodCount, 5);
for i = 1:interictalPeriodCount
    interictal{i} = interictalPeriod(epileptiformEventTimes(i,2)*frequency:epileptiformEventTimes(i+1,1)*frequency);    %Make vectors based on adjusted times to errors made by detection algorithm
    interictal{i} (interictal{i} == -1) = [];   %remove any spikes, artifacts or light pulses, that were marked for removal, during the interictal period 
    if length(interictal{i}) < (minInterictalPeriod*frequency)
        interictal{i} = -1; %This is a marker to ignore the interictal period below the minimum; I only want to analyze periods larger than 10 s
    end
    %Characterize baseline features from absolute value of the filtered/processed data
    interictal{i,4} = mean(interictal{i}); %Average
    interictal{i,5} = std(interictal{i}); %Standard Deviation
%     figure
%     plot (interictal{i})
%     title(sprintf('interictal period #%d. Sigma:%.4f', i, interictal{i,3}))
end

%Locate and delete the interictal period less than minimum Interictal Period secs
indexDelete = find ([interictal{:,4}] == -1); %locate 
interictal(indexDelete,:)=[]; %Delete
clear indexDelete


%Creating powerpoint slide
if figureInput == 1
    
isOpen  = exportToPPTX();
if ~isempty(isOpen)
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('new','Dimensions',[12 6], ...
    'Title','Epileptiform Event Detector V8.5', ...
    'Author','Michael Chang', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');
%Add New Slide
exportToPPTX('addslide');
exportToPPTX('addtext', 'Frequency Context of Epileptiform Events detected', 'Position',[2 1 8 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 36);
exportToPPTX('addtext', sprintf('File: %s', FileName), 'Position',[3 3 6 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
exportToPPTX('addtext', 'By: Michael Chang, Liam Long, and Kramay Patel', 'Position',[4 4 4 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
%Add New Slide
exportToPPTX('addslide');
exportToPPTX('addtext', 'Legend', 'Position',[0 0 4 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 24);
text = 'Note: If the dominant frequency detected was >50, it was ignored and replaced with 0 Hz (Bandpass filter set 1-50 Hz)';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 1 6 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = sprintf('Nyquist frequency (%.3f Hz) is the maximum valid frequency (sampling frequency/2)', frequency/2);
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 2 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = 'Rayleigh frequency: 1/windowSize (Hz), is the minimum frequency that can be resolved from signal';
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 3 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = sprintf('The window size used is %.1f s. so the minimum valid frequency is %.1f Hz', windowSize, 1/windowSize);
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 4 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
text = sprintf('Data was bandpass filtered (%s). Accordingly, if the dominant frequency detected was above %.0f Hz, it was considered invalid and reported as 0 Hz', filter, fc);
exportToPPTX('addtext', sprintf('%s',text), 'Position',[0 5 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 16);


% %Part D: Analysis
% %Locate the interictal with the lowest sigma, use as baseline
% [~, indexMin] = min ([interictal{:,5}]); %locate 
% 
% %Plot for your records
% i = indexMin;
% figHandle = figure;
% set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
% set(gcf, 'Position', get(0, 'Screensize'));
% subplot (2,1,1)
% plot (interictal{i})
% title(sprintf('Interictal period with lowest Sigma selected to be Baseline | Interictal Period #%d. Sigma: %.4f ', i, interictal{i,3}))
% ylabel('Voltage Activity (mV)')
% xlabel(sprintf('data points (Sampling Rate: %d Hz)', frequency))
% subplot (2,1,2)
% histogram(interictal{i})
% title(sprintf('Distribution of voltage activity from Interictal Period #%d', i))
% ylabel('Count (Frequency)')
% xlabel('Size of Voltage Activity (mV)')
% 
% %Export figures to .pptx
% exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
% exportToPPTX('addpicture',figHandle);
% close(figHandle)

end


%Calculate Frequency Content of Epileptiform Events
indexEvents = find(events(:,3) > windowSize);
p = round(frequency);   %to calculate the average freq/sec 
for i = indexEvents'
    %Event Vector
    eventVector = epileptiformEvent{i, 1};
    
    %Time Vector
    timeVector = (0:(length(eventVector)- 1))/frequency;
    timeVector = timeVector';    
    
    %Make filterbank | Focus on a specific frequency range 
    fb = cwtfilterbank('SignalLength', numel(eventVector), 'SamplingFrequency', frequency, 'FrequencyLimits', [0 fc], 'Wavelet','amor');
    %Perform continuous wavelet transform to calculat frequency content
    [wt, f] = cwt(eventVector, 'FilterBank',fb);
    p = abs(wt);    %Calculate power
    
    %Dominant Frequency at each time point 
    [maxS, idx] = max(p);        
    maxFreq = f(idx);   %finding the frequency with the maximum PSD
        
    
    %Calculate the average frequency per second
%     p = round(frequency);
    x = timeVector;
    x2 = maxFreq;

    S = numel(x);
    xx = reshape(x(1:S-mod(S,p)),p,[]);
    xx2 = reshape(x2(1:S-mod(S,p)),p,[]);
    
    y = sum(xx,1).'/p;
    y(:,2) = sum(xx2,1).'/p;
    
    %store the averaged max frequency content of each event 
    epileptiformEvent{i, 2} = y;
    epileptiformEvent{i, 3} = y(:,2);   
    
    clear x x2 S xx xx2 y
    
    %decipher
    [label, classification] = decipher (events,i);
    
    if figureInput == 1
   
    %Plot Figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('%s Event #%d (shading flat)', label, i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));

    s1=subplot (3,1,1);
    plot (timeVector, eventVector)  %Plot ictal event
        hold on
    plot (timeVector(round(windowSize*frequency)), eventVector(round(windowSize*frequency)), 'ro', 'color', 'black', 'MarkerFaceColor', 'green')    %SLE onset
    plot (timeVector(round(numel(eventVector)-(windowSize*frequency))), eventVector(round(numel(eventVector)-(windowSize*frequency))), 'ro', 'color', 'black', 'MarkerFaceColor', 'red')    %SLE offset
    title (sprintf('LFP Bandpass Filtered (%s), %s Event #%d   |   Treatment Group:%d', filter, label, i, events(i,4)))
    xlabel('Time (sec)')
    ylabel('Voltage (mV)')
    axis tight   
    
    s2=subplot (3,1,2);
    contourf(timeVector, f,p, 'edgecolor', 'none')
    shading flat
    c = colorbar;
    c.Label.String = 'Power Spectral Density (dB/Hz)';  %Originally I wrote 'Power (dB)', but I think I've been calculating the PSD
    title (sprintf('Frequency Content (PSD) of %s Event #%d. Michaels Algorithm detected: %s', label, i, classification))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')         
    
    %reposition scale bar
    s1Pos = get(s1,'position');
    s2Pos = get(s2,'position');
    s2Pos(3:4) = [s1Pos(3:4)];
    set(s2,'position',s2Pos);
    
    subplot (3,1,3)
    plot(timeVector,maxFreq) 
    title (sprintf('Dominant Frequency over duration of %s Event #%d', label, i))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    axis tight
    
     %Export figures to .pptx
     exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
     exportToPPTX('addpicture',figHandle);
     close(figHandle)
     
    end    
end
% 
% if figureInput == 1
% %Add New Slide
% exportToPPTX('addslide');
% exportToPPTX('addtext', 'Interictal periods used as Baseline Segments to normalized the frequency', 'Position',[2 1 8 2],...
%              'Horiz','center', 'Vert','middle', 'FontSize', 36);
% exportToPPTX('addtext', sprintf('File: %s', FileName), 'Position',[3 3 6 2],...
%              'Horiz','center', 'Vert','middle', 'FontSize', 20);
% text = 'The spectrogram of all the baselines, only the center portion equal to the window size was used to calculate PSD';
% exportToPPTX('addtext', sprintf('%s', text), 'Position',[4 4 4 2],...
%              'Horiz','center', 'Vert','middle', 'FontSize', 20);
% end
% 
% %Calculate Frequency Content of Baseline (interictal period)
% [nr, ~] = size (interictal);   %Count how many interictal periods there are, "nr"
% 
% for i = 1:nr
%     %Event Vector
%     eventVector = interictal{i, 1};
%     
%     %Time Vector
%     timeVector = (0:(length(eventVector)- 1))/frequency;
%     timeVector = timeVector';    
%     
%     %Frequency content of baseline event 
%     [s,f,t,p] = spectrogram (eventVector, round(windowSize*frequency), round(windowOverlap*frequency), 2.^nextpow2(windowSize*frequency), frequency, 'yaxis');
%     
%     %Dominant Frequency at each time point | NEW
%     [maxS, idx] = max(p);    
%     maxFreq = f(idx);
%     indexInvalid = find (maxFreq > fc);  %Ignore all frequency
%     maxFreq(indexInvalid) = 0;
%     
%     %store the max frequency content of each event 
%     interictal{i, 2} = t;
%     interictal{i, 3} = maxFreq;    
%     
%     if figureInput == 1
%     %decipher
%     label = 'baseline';
%     classification = 'baseline';   
%     
%     %Plot Figures
%     figHandle = figure;
%     set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('%s Event #%d', label, i)); %select the name you want
%     set(gcf, 'Position', get(0, 'Screensize'));
% 
%     subplot (3,1,1)
%     plot (timeVector(round(windowOverlap*frequency):round(t(end)*frequency)), eventVector(round(windowOverlap*frequency):round(t(end)*frequency)))
%     hold on
%     plot (timeVector(round(windowSize*frequency)), eventVector(round(windowSize*frequency)), 'ro', 'color', 'black', 'MarkerFaceColor', 'green')    %SLE onset
%     plot (timeVector(round(numel(eventVector)-(windowSize*frequency))), eventVector(round(numel(eventVector)-(windowSize*frequency))), 'ro', 'color', 'black', 'MarkerFaceColor', 'red')    %SLE offset
%     title (sprintf('LFP Bandpass Filtered (%s), %s Event #%d', filter, label, i))
%     xlabel('Time (sec)')
%     ylabel('Voltage (mV)')
%     axis tight   
%           
%     subplot (3,1,2)
%     imagesc(t,f,10*log10(p))
%     c = colorbar;
%     c.Label.String = 'Power (dB)';  
%     ylim([0 100])
%     title (sprintf('Frequency Content (PSD) of %s Event #%d. Michaels Algorithm detected: %s', label, i, classification))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (sec)')   
%     
%     subplot (3,1,3)
%     plot(t,maxFreq) 
%     title (sprintf('Dominant Frequency over duration of %s Event #%d', label, i))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (sec)')
%     axis tight   
%     ylim  ([0 70])            
%     
%     %Export figures to .pptx
%     exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
%     exportToPPTX('addpicture',figHandle);
%     close(figHandle)
%     end    
% end


%% save and close the .PPTX
if figureInput == 1
% % subtitle = guiInput{5};
% excelFileName = FileName(1:end-4);
% exportToPPTX('saveandclose',sprintf('%s(%s)', excelFileName, subtitle));
exportToPPTX('saveandclose',sprintf('%s', FileName));
end

%Closing Message
fprintf(1,'\nFrequency Content Analysis is Complete.\n')


   