function [indexEvents,widestGapThreshold ] = findThresholdSLE(featureSet,plotGraph)
%[indexEvents,widestGapThreshold ] = findThresholdSLE(featureSet,plotGraph)
%   This function classifies a populations of events into two distinct
%   populations based one feature set. This function is based on k-means
%   clustering and the widest gap method. The function will sort the
%   feature set and seperate the population based on where the widest gap
%   (in the value of the feature set) is observed, at or below the gap
%   created by k-means clustering. The populations of events with a larger
%   value (of the feature) is classied as SLEs (ictal events), and the
%   populatin with a lower value is classied as IIEs (interictal events). 
%   Algorithm is based on k-means clustering. K-means clustering is run 25x
%   on the feature set and uses the lowest threshold found. The
%   algorithm then looks for the widest gap below the lowest threshold and
%   uses that as the threshold that segments IIEs from SLEs. This function
%   find the most liberal threshold possible. There is a feature to plot
%   a graph of the feature sets and the thresholds found by the algorithm. 
%   Set the second input variable to '1' to plot graph for trouble
%   shooting.


if nargin == 1  
    plotGraph = 0;  %1 = yes; 0 = no
end   

%set variables
k = 2;  %number of cluster (IIE and SLE), set a apriori
featureSet = sort(featureSet);

%% Find lowest threshold
for i = 1:25
    %k-means clustering
%     idx = kmeans(featureSet, k, 'Start', 'cluster');      
    idx = kmeans(featureSet, k);      
    %find threshold
    group1 = featureSet(idx ==1); 
    group2 = featureSet(idx ==2); 
    if mean(group1)>mean(group2)
        upperLimit = min(group1);  %upper limit of threshold     
        lowerLimit = max(group2);  %lower limit of threshold
    else
        upperLimit = min(group2);  %upper limit of threshold     
        lowerLimit = max(group1);  %lower limit of threshold
    end    
    algoThreshold = (upperLimit + lowerLimit)/2; %assumed to be in between    
    %store threshold found
    thresholdStorage (i) = algoThreshold;            
end

%take lowest threshold calculated
minAlgoThreshold = min(thresholdStorage);

%find index of upper limit
indexUpperLimit = find(featureSet>minAlgoThreshold, 1); %the first index larger than the threshold - this is the upper limit
%% Find widest gap, (south of the algo determined threshold)

gap = diff(featureSet(:,1));    %find the size of gaps
[B,I] = sort (gap, 'ascend');   %sort them in order of smallest to largest

for i = numel(gap):-1:1         %find the widest gap below the algo-determined threshold
    if I(i) <= indexUpperLimit
        break
    end
end

indexWidestGap = I(i); 

indexWidestGapLowerLimit = indexWidestGap; 
indexWidestGapUpperLimit = indexWidestGapLowerLimit + 1; 

widestGapLowerLimit = featureSet(indexWidestGapLowerLimit); 
widestGapUpperLimit = featureSet(indexWidestGapUpperLimit); 

widestGapThreshold = (widestGapUpperLimit + widestGapLowerLimit)/2;

%segment events into SLE and IIE, using widest gap threshold
indexEvents = featureSet>widestGapThreshold;

%% plot figure

if plotGraph == 1    
    figure;
    gscatter(featureSet , featureSet, indexEvents);    %plot scatter plot
    hold on
    plot ([minAlgoThreshold minAlgoThreshold], ylim); %plot vertical line at 1    
    plot ([widestGapThreshold widestGapThreshold], ylim); %plot vertical line at 1
    legend ('IIEs', 'SLEs', 'Algo Threshold', 'widest-Gap Threshold')
    title ('Unsupervised SLE threshold for classification')
end

end

