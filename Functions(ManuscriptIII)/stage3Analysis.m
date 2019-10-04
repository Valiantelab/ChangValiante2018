function [resultsDuration, FigA] = stage3Analysis (durationControl,type, figure)
%Stage3Analysis.m is a function designed to perform the series of
%statistical analysis (paarametric or non-parametric) on an array of 
%numbers and determine if it is normally distributed based on the 
%Anderson-Darling Test. A figure of the histrogram can optionally plotted
%Output: an array of the results (column 1: mean/median; column 2: std/IQR;
%column 3: p-value from AD Test; If dealing with non-parametric data, 
%column 4 is the 2nd quartile and column 5 is the 4th quartile.
%Input: durationControl is an array of values, type is the type of data 
%'parametric' or 'nonparametric', and figure, 'figure' if you want a 
%figure, leave blank if otherwise. 


%If default values are not provided
if nargin <2
    type = 'parametric';
    figure = 'no figure';
end   

if nargin <3
    figure = 'no figure';
end   

%Statistical summary
switch type 
    case 'parametric'
        resultsDuration(1,1) = mean(durationControl);
        resultsDuration(1,2) = std(durationControl);
    case 'nonparametric'           
        resultsDuration(1,1) = median(durationControl); %Median value
        resultsDuration(1,4) = prctile(durationControl,25);  %25th percentile
        resultsDuration(1,5) = prctile(durationControl,75);  %75th percentile
        resultsDuration(1,2) = resultsDuration(1,5) - resultsDuration(1,4); %Interquartile Range (IQR)
    otherwise
        fprintf(2,'\nSecond input for stage2Analysis.m should be "parametric" or "nonparametric"\n')    %error meesage
end

%AD Test for Normality
if numel(durationControl)>3 %Need >3 samples to perform AD test for normality    
    [h,p] = adtest(durationControl); 
    if h == 0
        msg = 'normally distributed';
    else
        msg = 'not normally distributed';
    end
    resultsDuration(1,3) = p;
else
    resultsDuration(1,3) = NaN;
end

%Plot figure (optional)
switch figure 
    case 'figure'
        %plot distribution for visualization
        FigA = histfit(durationControl, 20);
        title(sprintf('Feature of ictal events are %s, p = %.2f', msg, p))
        xlabel ('Feature')
        ylabel ('Frequency (# of ictal events)')
    otherwise        
end

end
