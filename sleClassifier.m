function [indexSLE,featureThreshold ] = sleClassifier(featureSet,plotGraph)
%sleClassifer will be used to classify epileptiforms events as SLEs (ictus)
%   Algorithm is k-means clustering

if nargin == 1  
    plotGraph = 0;
end   

%set variables
k = 2;  %number of cluster, set a apriori

idx = kmeans(featureSet, k, 'Start', 'cluster');

% finding the threshold
groupSLE = featureSet(idx ==1); %SLEs are labelled 1
groupIIE = featureSet(idx ==2); %SLEs are labelled 2
upperLimit = min(groupSLE);  %upper limit of threshold
lowerLimit = max(groupIIE);  %lower limit of threshold
featureThreshold = (upperLimit + lowerLimit)/2; %assumed to be in between

indexSLE = featureSet>featureThreshold;

if plotGraph == 1    
    figure;
    gscatter(featureSet , featureSet, idx);    %plot scatter plot
    hold on
    plot ([featureThreshold featureThreshold], ylim); %plot vertical line at 1
end

end

