%% Stage 3: Analyze the file
% Author: Michael Chang
% Run this file after Stage 2 to analyze the results and do
% additional analysis to the detected events. This creats the time vector,
% LFP time series, LED if there is light, and filters the data using a
% bandpass filter (1-50 Hz) and a low pass filter (@68 Hz)


%Add all subfolders in working directory to the path.
addpath(genpath(pwd));  

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
histfit(durationControl, 20);
title(sprintf('Duration of ictal events are %s, p = %.2f', msg, p))
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsDuration(1,1) = mean(durationControl);
resultsDuration(1,2) = std(durationControl);

%Test Condition
if numel(durationTest) > 3 %Need 4 or more samples to perform AD test for normality
    [h, p] = adtest(durationTest);
    if h == 0
        msg = 'normally distributed';
    else
        msg = 'not normally distributed';
    end
    resultsDuration(2,3) = p;
else
    resultsDuration(2,3) = NaN;
end
%plot distribution for visualization
histfit(durationTest, 20);
title(sprintf('Duration of ictal events are %s, p = %.2f', msg, p))
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsDuration(2,1) = mean (durationTest);
resultsDuration(2,2) = std(durationTest);

%Posttest Conditions
if numel(durationPosttest) >2
    if numel(durationPosttest) >3
        [h,p] = adtest(durationPosttest);
        if h == 0
            msg = 'normally distributed';
        else
            msg = 'not normally distributed';
        end
        resultsDuration(3,3) = p;
    else
        resultsDuration(3,3) = NaN;
    end
%plot distribution for visualization
histfit(durationPosttest, 20);
title(sprintf('Duration of ictal events are %s, p = %.2f', msg, p))
xlabel ('Duration (s)')
ylabel ('Frequency (# of ictal events)')
%Staistics
resultsDuration (3,1) = mean(durationPosttest);
resultsDuration (3,2) = std(durationPosttest);
end

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
p_value_KS_duration = p;

%one-way ANOVA
[p,tbl,stats] = anova1(durationMatrix);
title('Box Plot of ictal event duration from different time periods')
xlabel ('Time Period')
ylabel ('Duration of ictal events (s)')
tbl_1ANOVA_duration = tbl;
%multiple comparisons, Tukey-Kramer Method
c = multcompare(stats);
c_duration = c;
%Kruskal Wallis
% p = kruskalwallis(durationMatrix);
% title('Boxplot: duration of ictal events from different treatment groups')
% xlabel ('Treatment Group')
% ylabel ('Duration (s)')


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
if numel(intensityTest)>3   %Need 4 or more samples to perform AD test for normality
    [h, p] = adtest(intensityTest);
    resultsIntensity(2,3) = p;
else
    resultsIntensity(2,3) = NaN;
end
%Plot distribution for visualization
histfit(intensityTest, 20);
title('Frequency histogram for intensity of ictal event')
xlabel ('Intensity (mW^2/s)')
ylabel ('Frequency (# of ictal events)')
%Statistics
resultsIntensity(2,1) = mean(intensityTest);
resultsIntensity(2,2) = std(intensityTest);

%Posttest Conditions
if numel(intensityPosttest)>2
    if numel(intensityPosttest)>3
        [h,p] = adtest(intensityPosttest);
        resultsIntensity(3,3) = p;
    else 
        resultsIntensity(3,3) = NaN;
    end    
    %Plot distribution for visualization
    histfit(intensityPosttest, 20);
    title('Frequency histogram for intensity of ictal event')
    xlabel ('Intensity (mW^2/s)')
    ylabel ('Frequency (# of ictal events)')
    %Statistics
    resultsIntensity(3,1) = mean(intensityPosttest);
    resultsIntensity(3,2) = std(intensityPosttest);
end

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
p_value_KS_intensity = p;


%one-way ANOVA
[p,tbl,stats] = anova1(intensityMatrix);
title('Box Plot of ictal event intensity from different time periods')
xlabel ('Treatment Condition')
ylabel ('intensity of ictal events (mV^2/s)')
tbl_1ANOVA_intensity = tbl;
%Multiple Comparisons, Tukey-Kramer Method
c = multcompare(stats);
c_intensity = c;
% %Kruskal Wallis
% p = kruskalwallis(intensityMatrix);
% title('Boxplot intensity of ictal events from different treatment groups')
% xlabel ('Treatment Group')
% ylabel ('intensity (e-5)')

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

%Figures for Visual Analysis
    FigE=figure;
    set(gcf,'Name','Control', 'NumberTitle', 'off');
    circ_plot(thetaControl,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Control Condition, p = %.3f', resultsTheta(1,1)));

    FigF=figure;
    set(gcf,'Name','Test','NumberTitle', 'off');
    circ_plot(thetaTest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Test Condition, p = %.3f', resultsTheta(2,1)));

if numel(thetaPosttest)>2
    resultsTheta(3,1)=circ_vtest(thetaPosttest,0);
    %Figures for visual analysis
    FigG=figure;
    set(gcf,'Name','Post-Test','NumberTitle', 'off');
    circ_plot(thetaPosttest,'hist',[],50,false,true,'linewidth',2,'color','r');
    title (sprintf('Post-Test Condition, p = %.3f',resultsTheta(3,1)));
end

%Combine all the results
result = horzcat(resultsDuration,resultsIntensity,resultsTheta);
%Label treatment groups
treatmentGroups = [1:3]';
%ictal events # in each group
n(1,1)=numel(thetaControl);
n(2,1)=numel(thetaTest);
n(3,1)=numel(thetaPosttest);

%% Write results to .xls 

excelFileName = 'result_acidosis.xlsx';
sheetName = FileName(1:8);

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

II = 'n';

J = 'KS Test, 1 vs 2';
K = 'one-way ANOVA, duration';
M = 'Multiple Comparison (Tukey-Kramer method), duration';
N = 'one-way ANOVA, intensity';
O = 'Multiple Comparison (Tukey-Kramer method), intensity';

P = 'Group';
Q = 'p-value';

%Write General Results
    subtitle1 = {A, B, C, D, E, F, G, H, I, II};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A1');
    xlswrite(sprintf('%s',excelFileName),treatmentGroups,sprintf('%s',sheetName),'A2');
    xlswrite(sprintf('%s',excelFileName),result,sprintf('%s',sheetName),'B2');
    xlswrite(sprintf('%s',excelFileName),n,sprintf('%s',sheetName),'J2');
%Write KS Test results
    subtitle1 = {J};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A6');
    xlswrite(sprintf('%s',excelFileName),p_value_KS_duration,sprintf('%s',sheetName),'B6');
    xlswrite(sprintf('%s',excelFileName),p_value_KS_intensity,sprintf('%s',sheetName),'E6');

% if numel(durationPosttest)>2
%Write one-way ANOVA results, duration
    subtitle1 = {K};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A8');
    xlswrite(sprintf('%s',excelFileName),tbl_1ANOVA_duration,sprintf('%s',sheetName),'A9');
%Write multiple comparison (Tukey-Kramer), duration
    subtitle1 = {M};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A14');
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'A15'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'B15'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{Q},sprintf('%s',sheetName),'F15'); %p-value subtitle
    xlswrite(sprintf('%s',excelFileName),c_duration,sprintf('%s',sheetName),'A16');
%Write one-way ANOVA results, intensity
    subtitle1 = {N};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A20');
    xlswrite(sprintf('%s',excelFileName),tbl_1ANOVA_intensity,sprintf('%s',sheetName),'A21');
%Write multiple comparison (Tukey-Kramer), duration
    subtitle1 = {O};
    xlswrite(sprintf('%s',excelFileName),subtitle1,sprintf('%s',sheetName),'A26');
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'A27'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{P},sprintf('%s',sheetName),'B27'); %Group subtitle
    xlswrite(sprintf('%s',excelFileName),{Q},sprintf('%s',sheetName),'F27'); %p-value subtitle
    xlswrite(sprintf('%s',excelFileName),c_intensity,sprintf('%s',sheetName),'A28');
% end

    