%% Stage 3: Analyze the file
% Author: Michael Chang
% Run this file after Stage 2 to analyze the results and do
% additional analysis to the detected events. This creats the time vector,
% LFP time series, LED if there is light, and filters the data using a
% bandpass filter (1-50 Hz) and a low pass filter (@68 Hz)

%% Write results to .xls 

excelFileName = 'result.xlsx';
sheetName = FileName(1:8);
treatmentGroups = [1:3]';

%set subtitle
A = 'Treatment Group';
B = 'Duration (s), average';
C = 'Duration (s), Std';
D = 'AD test, normality';
E = 'Intensity (mV^2/s), average';
F = 'Intensity (mV^2/s), std';
G = 'AD test, normality';
H = 'Light-triggered';
I = 'Dominant Frequency';

J = 'KS Test, 1 vs 2';
K = 'one-way ANOVA, duration';
M = 'Multiple Comparison (Tukey-Kramer method), duration';
N = 'one-way ANOVA, intensity';
O = 'Multiple Comparison (Tukey-Kramer method), intensity';

P = 'Group';
Q = 'p-value';

    subtitle1 = {A, B, C, D, E, F, G, H, I};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A1');
    xlswrite(sprintf('%s',excelFileName),treatmentGroups,sprintf('%s',sheetName),'A2');
    xlswrite(sprintf('%s',excelFileName),result,sprintf('%s',sheetName),'B2');

%% Duration
%Organize groups 
feature = 3;    %column #, i.e., duration is the 3rd column
durationControl = events(events(:,4)==1,feature);
durationTest = events(events(:,4)==2,feature);
durationPosttest = events(events(:,4)==3,feature);

%Control Condition
[h,p] = adtest(durationControl); %AD test for normality
if h == 0
    msg = 'normally distributed';
else
    msg = 'not normally distributed';
end
resultsDuration(1,3) = p;
%plot distribution for visualization
histfit(durationControl, 20)
title(sprintf('Duration of ictal events are %s, p = %.2f', msg, p))
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsDuration(1,1) = mean(durationControl);
resultsDuration(1,2) = std(durationControl);

%Test Condition
[h, p] = adtest(durationTest)
if h == 0
    msg = 'normally distributed';
else
    msg = 'not normally distributed';
end
resultsDuration(2,3) = p;
%plot distribution for visualization
histfit(durationTest, 20)
title(sprintf('Duration of ictal events are %s, p = %.2f', msg, p))
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsDuration(2,1) = mean (durationTest);
resultsDuration(2,2) = std(durationTest);

%Posttest Conditions
[h,p] = adtest(durationPosttest)
if h == 0
    msg = 'normally distributed';
else
    msg = 'not normally distributed';
end
resultsDuration(3,3) = p;
%plot distribution for visualization
histfit(durationPosttest, 20)
title(sprintf('Duration of ictal events are %s, p = %.2f', msg, p))
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Staistics
resultsDuration (3,1) = mean(durationPosttest);
resultsDuration (3,2) = std(durationPosttest);

%duration Matrix
durationMatrix(1:numel(durationControl),1) = durationControl;
durationMatrix(1:numel(durationTest),2) = durationTest;
durationMatrix(1:numel(durationPosttest),3) = durationPosttest;
durationMatrix(durationMatrix==0) = NaN;

%Analysis, comparison
%Independent sample Student's T-test
[h,p,ci,stats] = ttest2(durationControl, durationTest);
%2-sample KS Test
[h,p] = kstest2(durationControl, durationTest);
%one-way ANOVA
[p,tbl,stats] = anova1(durationMatrix);
title('Box Plot of ictal event duration from different time periods')
xlabel ('Time Period')
ylabel ('Duration of ictal events (s)')
%multiple comparisons, Tukey-Kramer Method
c = multcompare(stats);
%Kruskal Wallis
p = kruskalwallis(durationMatrix);
title('Boxplot: duration of ictal events from different treatment groups')
xlabel ('Treatment Group')
ylabel ('Duration (s)')

%% Intensity
%Organize groups
feature = 5;
intensityControl = events(events(:,4)==1,feature);
intensityTest = events(events(:,4)==2,feature);
intensityPosttest = events(events(:,4)==3,feature);

%Control Condition
[h,p] = adtest(intensityControl);
if h == 0
    msg = 'normally distributed';
else
    msg = 'not normally distributed';
end
resultsIntensity(1,3) = p;
%plot distribution for visualization
histfit(intensityControl, 20);
title(sprintf('Intensity of ictal events are %s, p = %.2f', msg, p))
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsIntensity(1,1) = mean(intensityControl);
resultsIntensity(1,2) = std(intensityControl);


%Test Condition
[h, p] = adtest(intensityTest);
resultsIntensity(2,3) = p;
%Plot distribution for visualization
histfit(intensityTest, 20);
title('Frequency histogram for intensity of ictal event')
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsIntensity(2,1) = mean(intensityTest);
resultsIntensity(2,2) = std(intensityTest);

%Posttest Conditions
[h,p] = adtest(intensityPosttest);
resultsIntensity(3,3) = p;
%Plot distribution for visualization
histfit(intensityPosttest, 20);
title('Frequency histogram for intensity of ictal event')
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsIntensity(3,1) = mean(intensityPosttest);
resultsIntensity(3,2) = std(intensityPosttest);

%intensity Matrix
intensityMatrix(1:numel(intensityControl),1) = intensityControl;
intensityMatrix(1:numel(intensityTest),2) = intensityTest;
intensityMatrix(1:numel(intensityPosttest),3) = intensityPosttest;
intensityMatrix(intensityMatrix==0) = NaN;

%Analysis, comparison
%Independent sample Student's T-test
[h,p,ci,stats] = ttest2(intensityControl, intensityTest);
%2-sample KS Test
[h,p] = kstest2 (intensityControl, intensityTest);
%one-way ANOVA
[p,tbl,stats] = anova1(intensityMatrix)
title('Box Plot of ictal event intensity from different time periods')
xlabel ('Treatment Condition')
ylabel ('intensity of ictal events (mV^2/s)')
%Multiple Comparisons, Tukey-Kramer Method
c = multcompare(stats);
%Kruskal Wallis
p = kruskalwallis(intensityMatrix);
title('Boxplot intensity of ictal events from different treatment groups')
xlabel ('Treatment Group')
ylabel ('intensity (e-5)')


%% Circular Variance and Plots of Ictal Event with photosimulation    
%Calculate Theta
for i = 1:numel(events(:,1))
    events(i, 29) = events(i, 26)/events(i, 27) * (2*pi);
end
%Organize Group
feature = 29;
thetaControl = events(events(:,4)==1,feature);
thetaTest = events(events(:,4)==2,feature);
thetaPosttest = events(events(:,4)==3,feature);

% thetaControl=SLE(controlStart<SLE(:,1) & SLE(:,1)<controlEnd,feature);
% thetaTest=SLE(testStart<SLE(:,1) & SLE(:,1)<testEnd,feature);
% thetaPosttest=SLE(posttestStart<SLE(:,1) & SLE(:,1)<posttestEnd,feature);

%Analysis, light correlation?
resultsTheta(1,1)=circ_vtest(thetaControl,0);
resultsTheta(2,1)=circ_vtest(thetaTest,0);
resultsTheta(3,1)=circ_vtest(thetaPosttest,0);

%Figures for Visual Analysis
    FigE=figure;
    set(gcf,'Name','G. theta 4-Aminopyrimidine', 'NumberTitle', 'off');
    circ_plot(thetaControl,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Control Condition, p = %.3f', resultsTheta(1,1)));

    FigF=figure;
    set(gcf,'Name','I. theta Hepes-Buffered ACSF','NumberTitle', 'off');
    circ_plot(thetaTest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Test Condition, p = %.3f', resultsTheta(2,1)));
    
    FigG=figure;
    set(gcf,'Name','I. theta Hepes-Buffered ACSF','NumberTitle', 'off');
    circ_plot(thetaPosttest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Post-Test Condition, p = %.3f',resultsTheta(3,1)));

%Combine all the results
result = horzcat(resultsDuration,resultsIntensity,resultsTheta);