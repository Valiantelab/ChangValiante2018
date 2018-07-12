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
fs = 1000000/si; %Hz si is the sampling interval in microseconds from the metadata
t = (0: (length(x)-1))/fs;

%% Seperate .abf file signals into independent signals
LFP = x(:,1); 
LED = x(:,2);  %%To be used if you need to collect LED data in the future' switch column 1 to 2

%% import SLE onset and offset times from excel
%import excel times to create vectors
SLE_onset_raw= xlsread(excel_filename,1,'C2:C50');
SLE_offset_raw= xlsread(excel_filename,1,'D2:D50');
duration_raw = xlsread(excel_filename,1,'E2:E50');
onset_delay = xlsread(excel_filename,1,'B2:B50');
%Finding Baseline onset and offset times
baseline_onset_raw = SLE_onset_raw - duration_raw;  %Make sure there is baseline in front of 1st SLE
         
%% Convert seconds to corresponds with the original data points (LFP)
SLE_onset=round(SLE_onset_raw*10000);   
SLE_offset=round(SLE_offset_raw*10000);
baseline_onset=round(baseline_onset_raw*10000);


%% Creating test vector
%Preallocate cell arrays
SLE =cell(1, numel(SLE_onset));
SLETime =cell(1, numel(SLE_onset));

for i=1:numel(SLE_onset);
        clear y0_d x0_dd y1_d x1_dd control_distance_between_points ictal_event_distance_between_points n nb n2 nb2 light0_d light1_d

        %Ictal Event vector
        y1_d = LFP (SLE_onset(i):SLE_offset(i)); % y1_d is the decimated ictal event LFP vector                
        x1_dd = t (SLE_onset(i):SLE_offset(i)); % x1_dd is the decimated ictal event time vector
        
        %store the vectors                
        SLE{i} = y1_d;
        SLETime{i} = x1_dd;

end


%% Testing pwelch to find freqyenct content of seizure 
%(Taken from MatLab Website)

fsamp = 10000;
xx = SLE{2};    %datapoints for the 2nd SLE 
tt = SLETime{2};    %Time vector for the 2nd SLE
 
% Plot time-domain signal
subplot(2,1,1);
plot(tt, xx);
ylabel('Amplitude'); xlabel('Time (secs)');
axis tight;
title('Noisy Input Signal');

% Choose FFT size and calculate spectrum
Nfft = 1024;
[Pxx,f] = pwelch(xx,gausswin(Nfft),Nfft/2,Nfft,fsamp);

% Plot frequency spectrum
subplot(2,1,2);
plot(f,Pxx);
ylabel('PSD'); xlabel('Frequency (Hz)');
grid on;

% Get frequency estimate (spectral peak)
[~,loc] = max(Pxx);
FREQ_ESTIMATE = f(loc)
title(['Frequency estimate = ',num2str(FREQ_ESTIMATE),' Hz']);


'success; complete!'
