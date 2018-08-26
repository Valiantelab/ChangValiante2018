close all
clear all
clc

%% Test data to analyze
excel_filename = '13226009.xlsx'    % Enter the excel file name
    
%% Load 13226009(filtered).abf file                                            
[FileName,PathName] = uigetfile ('*.abf','pick 13226009(filtered).abf',...
    'F:');%Choose abf file
[x,si,metadata]=abfload ([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to
                                               %time value by dividing 10k
                                               
%% create time vector
fs = 1e6/si; %Hz si is the sampling interval in microseconds from the metadata
t = (0: (length(x)-1))/fs;

%% Seperate .abf file signals into independent signals
LFP = x(:,1); 
LED = x(:,2);  %%To be used if you need to collect LED data in the future' switch column 1 to 2


%% Testing pwelch to find freqyenct content of seizure 
%(Taken from MatLab Website)

i = 17;
k = 2;
SLE_Vector{i,k};
xx = data1(SLE_Vector{i,k});
time = t;
fs = frequency;

% Plot time-domain signal
figure(152); subplot(2,1,1);
plot(time, xx);
ylabel('Amplitude'); xlabel('Time (secs)');
axis tight;
title('Noisy Input Signal');

fftFunc = @(winData) pwelch_plot(winData, fs);

figure(101)
%third and fourth arguments specify window size and overlap respectively
%window size - length of each window (e.g. fs - > windows are 1 second long)
%overlap - how much consective windows overlap (e.g. fs/2 -> windows move 0.5s each fft)
[idx, results] = slideWindow(fftFunc, xx, fs,fs/2);

figure(102)
yyaxis left
plot(time, detrend(xx));
ylabel('Detrended LFP (mV)')
xlabel('Time(s)')
ylim([-2 2])
yyaxis right
scatter(time(idx),results);
ylabel('Highest Power Frequency (Hz)')

'success; complete!'

function [idx,result] = slideWindow(func, data, windowLength, overlap)

    if nargin == 3
        overlap = 0;
    end 
    
    windowEnd = windowLength;
    idx = [];
    result = [];
    while (windowEnd < length(data))
        idx = [idx  windowEnd];
        currentWindow = data(windowEnd-windowLength+1:windowEnd);
        
        result = [result func(currentWindow)];
        
        
        
        windowEnd = windowEnd + windowLength-overlap;
    end
    
    
    
end
    


% Choose FFT size and calculate spectrum
function FREQ_ESTIMATE = pwelch_plot(xx,fsamp)

    xx = detrend(xx);
    
    p = nextpow2(numel(xx));
    Nfft = 2^p;

    padding = Nfft - numel(xx);
    padding_zeros = zeros (padding, 1);

    xx = [xx; padding_zeros];   %xx is now padded with zeros


    %pwelch function
    [Pxx,f] = pwelch(xx);
    
    %plot time signal samples
    subplot(2,1,1);
    plot(xx);
    ylim([-0.5 0.5]);
    title('Time Signal Data');

    % Plot frequency spectrum
    subplot(2,1,2);
    plot (f*fsamp/2/pi,Pxx)
    %loglog(f,Pxx);
    ylabel('PSD'); xlabel('Frequency (Hz)');
    grid on;
    xlim ([0 100])

    % Get frequency estimate (spectral peak)
    [~,loc] = max(Pxx);
    FREQ_ESTIMATE = f(loc)*fsamp/2/pi;
    title(['Frequency estimate = ',num2str(FREQ_ESTIMATE),' Hz']);
    
    pause(0.01)
   
end

