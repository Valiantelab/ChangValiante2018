%% Stage 2: Process the File
% Author: Michael Chang
% Run this file after the detection algorithm to analyze the results and do
% additional analysis to the detected events. This creats the time vector,
% LFP time series, LED if there is light, and filters the data using a
% bandpass filter (1-50 Hz) and a low pass filter (@68 Hz)

%Create time vector
frequency = 1000000/samplingInterval; %Hz. si is the sampling interval in microseconds from the metadata
t = (0:(length(x)- 1))/frequency;
t = t';

%Seperate signals from .abf files
LFP = x(:,1);   %original LFP signal
if userInput(4)>0
    LED = x(:,userInput(4));   %light pulse signal, as defined by user's input via GUI
    onsetDelay = 0.13;  %seconds
    offsetDelay = 1.5;  %seconds
    lightpulse = LED > 1;
else
    LED =[];
    onsetDelay = [];
end

%Filter Bank
%Band Pass Filter
[b,a] = butter(2, ([1 50]/(frequency/2)), 'bandpass');  %Band pass filter
LFP_filteredBandPass = filtfilt (b,a,LFP);             %Bandpass filtered [1 - 50 Hz] singal; because of the 76 Hz noise above, also SLEs only have frequencies up to 20 Hz

%Low Pass Filter
fc = 68; % Cut off frequency
[b,a] = butter(4,fc/(frequency/2), 'low'); %Butterworth filter of order 4
LFP_filtered = filtfilt(b,a,LFP_filteredBandPass); %filtered signal

%Notch Filter
%@60 and @76 Hz