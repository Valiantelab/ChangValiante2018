close all
clear all
clc

%% Enter data to make comparison using the Two-sample Kolmogorov-Smirnov test to see if seizures changed power after drug application 

    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', 'C:\Users\User\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)');%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
                                       
    [excelFileName,excelPathName] = uigetfile('*.xlsx','Select the Excel file');
    sheetTitle = 'Experimental Protocol #1';
    [excelNum, excelText, excelRaw] = xlsread(strcat(excelPathName, excelFileName), sheetTitle);
    
    seizureOnsetTimes = excelNum(:,2);
    seizureOffsetTimes = excelNum(:,3);
    durations = excelNum(:,4);
    
    numElements = numel(find(~isnan(seizureOnsetTimes))); %number of actual numbers, excludes NaNs from import
    seizureOnsetTimes = seizureOnsetTimes(1:numElements);
    seizureOffsetTimes = seizureOffsetTimes(1:numElements);
    durations = durations(1:numElements);
                 
%% create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';
%% Seperate .abf file signals into independent signals
LFP = x(:,1); 
LED = x(:,2);  %%To be used if you need to collect LED data in the future' switch column 1 to 2

% [P] = pulse_seq(LED);
% light_pulse_onset_original = P.range(:,1);
% light_pulse_offset_original = P.range(:,2);


%% normalize the data
LFP_normalized = LFP - LFP(1);
%% plot graph of normalized  data 
figure;
set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
set(gcf,'Name','Overview of Data'); %select the name you want
set(gcf, 'Position', get(0, 'Screensize'));

lightpulse = LED > 1;

plot (t, LFP_normalized, 'k')
hold on
plot (t,lightpulse - 5)
for i = 1:length(seizureOnsetTimes)
    plot(t(int32(frequency*seizureOnsetTimes(i) + 1)), LFP_normalized(int32(frequency*seizureOnsetTimes(i) + 1)), 'r', 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 20)
    plot(t(int32(frequency*seizureOffsetTimes(i) + 1)), LFP_normalized(int32(frequency*seizureOffsetTimes(i) + 1)), 'c', 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 20)
end

title ('Overview of LFP recording (10000 points/s)');
ylabel ('LFP (mV)');
xlabel ('Time (s)');
%saveas(gcf,strcat(excelPathName, sprintf('%s Overview.png',FileName(1:length(FileName)-4))))
close

%% Line Lengths (of total ictal event for now) 
[L, LperSecond] = findLineLengths(frequency, seizureOnsetTimes, seizureOffsetTimes, t, LFP_normalized);

%% (Test) print values of variables of interest
% L
% LperSecond
output = [L LperSecond];
outputRange = sprintf('E2:F%d', numElements + 1);
xlswrite(strcat(excelPathName, excelFileName), output, sheetTitle, outputRange);
sprintf('%s Write Complete!', excelFileName)
