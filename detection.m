%Program: Epileptiform Activity Detector 
%Author: Michael Chang (michael.chang@live.ca), Fred Chen and Liam Long; 
%Copyright (c) 2018, Valiante Lab
%Version 3.0

%% Clear All
close all
clear all
clc

%% GUI to set thresholds
%Settings, request for user input on threshold
titleInput = 'Specify Detection Thresholds';
prompt1 = 'Epileptiform Spike Threshold: average + (6 x Sigma)';
prompt2 = 'Artifact Threshold: average + (100 x Sigma) ';
prompt3 = 'Figure: Yes (1) or No (0)';
prompt = {prompt1, prompt2, prompt3};
dims = [1 70];
definput = {'6', '100', '0'};
opts = 'on';
threshold_multiple = str2double(inputdlg(prompt,titleInput,dims,definput, opts));

%setting on distance between spikes, hard coded
distanceSpike = 1;  %distance between spikes (seconds)
distanceArtifact = 0.6; %distance between artifacts (seconds)

%% Load .abf and excel data
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', 'C:\Users\User\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)\1) abf files');%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
                                                                                         
%% create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%% Seperate signals from .abf files
LFP = x(:,1);   %original LFP signal
LED = x(:,2);   %light pulse signal

%% Data Processing 
%Center the LFP data
LFP_normalized = LFP - LFP(1);                                      %centered signal at 0, y-axis

%Lowpass butter filter [2Hz]
fc = 2; % Cut off frequency
[b,a] = butter(2,fc/(frequency/2)); % Butterworth filter of order 2
LFP_normalizedLowPassFiltered = filtfilt(b,a,LFP_normalized); % Will be the filtered signal

%Bandpass butter filter [1 - 100 Hz]
[b,a] = butter(2, [[1 100]/(frequency/2)], 'bandpass');
LFP_normalizedFiltered = filtfilt (b,a,LFP_normalized);             %Filtered signal

%Absolute value of the filtered data
AbsLFP_normalizedFiltered = abs(LFP_normalizedFiltered);            %1st derived signal

%Derivative of the filtered data (absolute value)
DiffLFP_normalizedFiltered = abs(diff(LFP_normalizedFiltered));     %2nd derived signal

% %Power of the derivative of the filtered data (absolute values)     
% powerFeature = (DiffLFP_normalizedFiltered).^2;                     %3rd derived signal

%% Detect potential events (epileptiform/artifacts) | Derivative Values
[epileptiformLocation, artifacts, locs_spike_1st] = detectEvents (DiffLFP_normalizedFiltered, frequency);

%remove potential events
for i = 1:size(epileptiformLocation,1)
AbsLFP_normalizedFiltered (epileptiformLocation (i,1):epileptiformLocation (i,2)) = [-1];
end

%remove artifacts
for i = 1:size(artifacts,1)
AbsLFP_normalizedFiltered (artifacts(i,1):artifacts(i,2)) = [-1];
end

%Isolate baseline recording
AbsLFP_normalizedFiltered (AbsLFP_normalizedFiltered == -1) = [];
AbsLFP_normalizedFilteredBaseline = AbsLFP_normalizedFiltered; %Rename

%Characterize baseline features from absolute value of the filtered data 
avgBaseline = mean(AbsLFP_normalizedFilteredBaseline); %Average
sigmaBaseline = std(AbsLFP_normalizedFilteredBaseline); %Standard Deviation

%% Detect events (epileptiform/artifacts) | Absolute Values

%Recreate the Absolute filtered LFP (1st derived signal) vector
AbsLFP_normalizedFiltered = abs(LFP_normalizedFiltered); %the LFP analyzed

%Define thresholds for detection, using inputs from GUI
minPeakHeight = avgBaseline+(threshold_multiple(1)*sigmaBaseline);      %threshold for epileptiform spike detection
minPeakDistance = distanceSpike*frequency;                              %minimum distance spikes must be apart
minArtifactHeight = avgBaseline+(threshold_multiple(2)*sigmaBaseline);  %threshold for artifact spike detection
minArtifactDistance = distanceArtifact*frequency;                       %minimum distance artifact spikes must be apart

%Detect events
[epileptiformLocation, artifacts, locs_spike_2nd] = detectEvents (AbsLFP_normalizedFiltered, frequency, minPeakHeight, minPeakDistance, minArtifactHeight, minArtifactDistance);

%% Finding event time 
epileptiformTime = [epileptiformLocation/frequency];

%% Classifier
SLE = epileptiformTime(epileptiformTime(:,3)>=10,:);
IIS = epileptiformTime(epileptiformTime(:,3)<10,:);

%% SLE: Determine exact onset and offset times | Power Feature
% Scan Low-Pass Filtered Power signal for precise onset/offset times
SLE_final = SLECrawler(LFP_normalizedFiltered, SLE, frequency, LED, 0.13, locs_spike_2nd);

%Store light-triggered events (s)
triggeredEvents = SLE_final(SLE_final(:,4)>0, 1);

%% Write to .xls
excelFileName = FileName(1:8);
A = 'Onset (s)';
B = 'Offset (s)';
C = 'Duration (s)';
D = 'Light-triggered (1 = yes)';

if isempty(SLE_final) == 0
    %Sheet 1 = SLE
    subtitle1 = {A, B, C, D};
    xlswrite(sprintf('%s(algo)',excelFileName),subtitle1,'SLE' ,'A1');
    xlswrite(sprintf('%s(algo)',excelFileName),SLE_final,'SLE' ,'A2');
else
    display ('No SLEs were detected.');
end
    
%Sheet 2 = IIS
subtitle2 = {A, B, C, D};
xlswrite(sprintf('%s(algo)',excelFileName),subtitle2,'IIS' ,'A1');
xlswrite(sprintf('%s(algo)',excelFileName),IIS,'IIS' ,'A2');

if isempty(artifacts) == 0
    %Sheet 3 = Artifacts
    subtitle3 = {A, B, C};
    xlswrite(sprintf('%s(algo)',excelFileName),subtitle3,'Artifacts','A1');
    xlswrite(sprintf('%s(algo)',excelFileName),artifacts/frequency,'Artifacts','A2');
else
    display ('No artifacts were detected.');
end


%% Optional: Plot Figures
if threshold_multiple(3) == 1   
    
%% Creating powerpoint slide
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end

exportToPPTX('new','Dimensions',[12 6], ...
    'Title','Epileptiform Detector V4.0', ...
    'Author','Michael Chang', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');

exportToPPTX('addslide');
exportToPPTX('addtext', 'SLE Events detected ', 'Position',[2 1 8 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 36);
exportToPPTX('addtext', sprintf('File: %s', FileName), 'Position',[3 3 6 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);
exportToPPTX('addtext', 'By: Michael Chang and Christopher Lucasius', 'Position',[4 4 4 2],...
             'Horiz','center', 'Vert','middle', 'FontSize', 20);     
         
exportToPPTX('addslide');
exportToPPTX('addtext', 'Legend', 'Position',[0 0 4 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 24);
exportToPPTX('addtext', 'Epileptiform spike is average + 6*SD of the baseline', 'Position',[0 1 6 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'Artifacts are average + 100*SD', 'Position',[0 2 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'SLE onset is the first peak in power (minimum 1/3 of the max amplitude spike)', 'Position',[0 3 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'SLE offset is when power returns below baseline/2', 'Position',[0 4 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 14);
exportToPPTX('addtext', 'Note: The event have only been shifted alone the y-axis to start at position 0', 'Position',[0 5 5 1],...
             'Horiz','left', 'Vert','middle', 'FontSize', 16);      


%% Plotting out detected SLEs with context | To figure out how off you are
data1 = LFP_normalized; %Time series to be plotted 

for i = 1:size(SLE,1)
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('V4.0 SLE #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
   
    time1 = single(SLE_final(i,1)*10000);
    time2 = single(SLE_final(i,2)*10000);
    rangeSLE = (time1:time2);
    if time1 >= 50001
        rangeOverview = (time1-50000:time2+50000);
    else
        rangeOverview = (1:time2+50000);
    end
        
    plot (t(rangeOverview),data1(rangeOverview))
    hold on
    plot (t(rangeSLE),data1(rangeSLE))     %SLE
    plot (t(time1), data1(time1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
    plot (t(time2), data1(time2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
    title (sprintf('LFP Recording, SLE #%d', i));
    ylabel ('mV');
    xlabel ('Time (sec)');
   
    exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
    exportToPPTX('addpicture',figHandle);      
    close(figHandle)
end
        
% save and close the .PPTX
newFile = exportToPPTX('saveandclose',sprintf(excelFileName)); 

%% plot entire recording 
figHandle = figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name', sprintf ('Overview of %s', FileName)); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

lightpulse = LED > 1;

subplot (3,1,1)
reduce_plot (t, LFP_normalized, 'k');
hold on
reduce_plot (t, lightpulse - 2);

%plot artifacts (red), found in 2nd search
for i = 1:numel(artifacts(:,1)) 
    reduce_plot (t(artifacts(i,1):artifacts(i,2)), LFP_normalized(artifacts(i,1):artifacts(i,2)), 'r');
end

%plot onset markers
for i=1:numel(epileptiformTime(:,1))
reduce_plot ((onsetTimes(i)), (LFP_normalized(epileptiformLocation(i))), 'o');
end

%plot offset markers
for i=1:numel(epileptiformTime(:,2))
reduce_plot ((offsetTimes(i)), (LFP_normalized(epileptiformLocation(i,2))), 'x');
end

title (sprintf ('Overview of LFP (10000 points/s), %s', FileName));
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,2) 
reduce_plot (t, AbsLFP_normalizedFiltered, 'b');
hold on
reduce_plot (t, lightpulse - 1);

%plot spikes (artifact removed)
for i=1:size(locs_spike_2nd,1)
plot (t(locs_spike_2nd(i,1)), (DiffLFP_normalizedFiltered(locs_spike_2nd(i,1))), 'x')
end

title ('Overview of filtered LFP (bandpass: 1 to 100 Hz)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');

subplot (3,1,3) 
reduce_plot (t(1:end-1), DiffLFP_normalizedFiltered, 'g');
hold on

% %plot onset markers
% for i=1:numel(locs_onset)
% plot (t(locs_spike(locs_onset(i))), (pks_spike(locs_onset(i))), 'o')
% end
 
%plot spikes 
for i=1:size(locs_spike_1st,1)
plot (t(locs_spike_1st(i,1)), (DiffLFP_normalizedFiltered(locs_spike_1st(i,1))), 'x')
end

% %plot artifacts, found in the 1st search
% for i=1:size(locs_artifact_1st(: ,1))
% plot (t(locs_artifact_1st(i,1)), (DiffLFP_normalizedFiltered(locs_artifact_1st(i,1))), 'o', 'MarkerSize', 12)
% end

% %plot offset markers
% for i=1:numel(locs_onset)
% plot ((offsetTimes(i)), (pks_spike(locs_onset(i))), 'x')
% end

title ('Peaks (o) in Derivative of filtered LFP');
ylabel ('Derivative (mV)');
xlabel ('Time (s)');

saveas(figHandle,sprintf('%s.png', FileName), 'png');

end

'successfully completed. Thank you for choosing to use The Epileptiform Detector'
