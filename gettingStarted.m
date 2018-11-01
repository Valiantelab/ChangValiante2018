%Load .abf file (raw data), analyze single file
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', inputdir);%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency)
    
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
[b,a] = butter(2, ([1 100]/(frequency/2)), 'bandpass');
LFP_filtered = filtfilt (b,a,LFP);             %Bandpass filtered [1 - 100 Hz] singal

%Plot Figure
figure;
plot(t, LFP_filtered)

