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

%% Time Lines for Comparisons

%Control Condition 
controlStart = 550;
controlEnd = 3270;
%Test Condition
testStart = 4975;
testEnd = 9300;
%Posttest Condition
posttestStart = 9360;
posttestEnd = 10140;

%% Duration
%Control Condition
feature = 3;
durationControl=SLE(controlStart<SLE(:,1) & SLE(:,1)<controlEnd,feature);
[h,p] = adtest(durationControl)
%plot distribution for visualization
histfit(durationControl, 20)
title('Frequency histogram for duration of ictal event')
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
median(durationControl);
prctile(durationControl, 25)
prctile(durationControl, 75)

%Test Condition
feature = 3;
durationTest=SLE(testStart<SLE(:,1) & SLE(:,1)<testEnd,feature);
[h, p] = adtest(durationTest)
%plot distribution for visualization
histfit(durationTest, 20)
title('Frequency histogram for duration of ictal event')
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
median(durationTest);
prctile(durationTest, 25)
prctile(durationTest, 75)

%Posttest Conditions
feature = 3;
durationPosttest=SLE(posttestStart<SLE(:,1) & SLE(:,1)<posttestEnd,feature);
[h,p] = adtest(durationPosttest)
%plot distribution for visualization
histfit(durationPosttest, 20)
title('Frequency histogram for duration of ictal event')
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Staistics
median(durationPosttest);
prctile(durationPosttest, 25)
prctile(durationPosttest, 75)

%duration Matrix
durationMatrix(1:numel(durationControl),1) = durationControl;
durationMatrix(1:numel(durationTest),2) = durationTest
durationMatrix(1:numel(durationPosttest),3) = durationPosttest
durationMatrix(durationMatrix==0) = NaN;

%Analysis, comparison
%Independent sample Student's T-test
[h,p,ci,stats] = ttest2(durationControl, durationTest)

%one-way ANOVA
[p,tbl,stats] = anova1(durationMatrix)
title('Box Plot of ictal event duration from different time periods')
xlabel ('Time Period')
ylabel ('Duration of ictal events (s)')

%Kruskal Wallis
p = kruskalwallis(durationMatrix)
title('Boxplot: duration of ictal events from different treatment groups')
xlabel ('Treatment Group')
ylabel ('Duration (s)')

%% Compare Intensity

%Control Condition
feature = 5;
intensityControl=SLE(controlStart<SLE(:,1) & SLE(:,1)<controlEnd,feature);
[h,p] = adtest(intensityControl);
%plot distribution for visualization
histfit(intensityControl, 20)
title('Frequency histogram for intensity of ictal event')
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
medianIntensity = median(intensityControl);
intensity25prctile = prctile(intensityControl, 25);
intensity75prctile = prctile(intensityControl, 75);


%Test Condition
feature = 5;
intensityTest=SLE(testStart<SLE(:,1) & SLE(:,1)<testEnd,feature);
[h, p] = adtest(intensityTest)
%Plot distribution for visualization
histfit(intensityTest, 20)
title('Frequency histogram for intensity of ictal event')
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
median(intensityTest)
prctile(intensityTest, 25)
prctile(intensityTest, 75)

%Posttest Conditions
feature = 5;
intensityPosttest=SLE(posttestStart<SLE(:,1) & SLE(:,1)<posttestEnd,feature);
[h,p] = adtest(intensityPosttest)
%Plot distribution for visualization
histfit(intensityPosttest, 20);
title('Frequency histogram for intensity of ictal event')
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
median(intensityPosttest);
prctile(intensityPosttest, 25)
prctile(intensityPosttest, 75)

%intensity Matrix
intensityMatrix(1:numel(intensityControl),1) = intensityControl;
intensityMatrix(1:numel(intensityTest),2) = intensityTest;
intensityMatrix(1:numel(intensityPosttest),3) = intensityPosttest;
intensityMatrix(intensityMatrix==0) = NaN;

%Analysis, comparison
%Independent sample Student's T-test
[h,p,ci,stats] = ttest2(intensityControl, intensityTest)
%one-way ANOVA
[p,tbl,stats] = anova1(intensityMatrix)
title('Box Plot of ictal event intensity from different time periods')
xlabel ('Time Period')
ylabel ('intensity of ictal events (mV^2/s)')


%Kruskal Wallis
p = kruskalwallis(intensityMatrix)
title('Boxplot intensity of ictal events from different treatment groups')
xlabel ('Treatment Group')
ylabel ('intensity (e-5)')


%% Circular Variance and Plots of Ictal Event with photosimulation
    
%Calculate Theta
for i = 1:numel(SLE(:,1))
    SLE(i, 17) = SLE(i, 15)/SLE(i, 16) * (2*pi);
end

feature = 17;
thetaControl=SLE(controlStart<SLE(:,1) & SLE(:,1)<controlEnd,feature);
thetaTest=SLE(testStart<SLE(:,1) & SLE(:,1)<testEnd,feature);
thetaPosttest=SLE(posttestStart<SLE(:,1) & SLE(:,1)<posttestEnd,feature);
  
  vtest_P_value_control=circ_vtest(thetaControl,0);
  vtest_P_value_test=circ_vtest(thetaTest,0);
  vtest_P_value_posttest=circ_vtest(thetaPosttest,0);
  
    FigE=figure;
    set(gcf,'Name','G. theta 4-Aminopyrimidine', 'NumberTitle', 'off');
    circ_plot(thetaControl,'hist',[],50,false,true,'linewidth',2,'color','r');
    title ('Control Condition');
    
    FigF=figure;
    set(gcf,'Name','I. theta Hepes-Buffered ACSF','NumberTitle', 'off');
    circ_plot(thetaTest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title ('Test Condition');
    
    FigG=figure;
    set(gcf,'Name','I. theta Hepes-Buffered ACSF','NumberTitle', 'off');
    circ_plot(thetaPosttest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title ('Post-Test Condition');

    