function [indexEvents,widestGapThreshold ] = thresholdFinder(featureSet,plotGraph)
%sleClassifer will be used to classify epileptiforms events as SLEs (ictus)
%   Algorithm is k-means clustering

if nargin == 1  
    plotGraph = 0;
end   

%set variables
k = 2;  %number of cluster, set a apriori
featureSet = sort(featureSet)
thresholdStorage = zeros(1,25);

%% Find lowest threshold
for i = 1:25
    %k-means clustering
    idx = kmeans(featureSet, k, 'Start', 'cluster');      
    %find threshold
    groupSLE = featureSet(idx ==1); %SLEs are labelled 1
    groupIIE = featureSet(idx ==2); %IIEs are labelled 2
    upperLimit = min(groupSLE);  %upper limit of threshold     
    lowerLimit = max(groupIIE);  %lower limit of threshold
    featureThreshold = (upperLimit + lowerLimit)/2; %assumed to be in between
    
%     indexUpperLimit = find(featureSet==upperLimit); %find index 
    %store threshold found
    thresholdStorage (i,1) = featureThreshold;    
    %store index
%     thresholdStorage (i,2) = indexUpperLimit;
    
end

%take lowest threshold calculated
minFeatureThreshold = min(thresholdStorage)

%find index of upper limit
indexUpperLimit = find(featureSet > minFeatureThreshold, 1); %the first index larger than the threshold - this is the upper limit
%% Find widest gap, (south of the algo determined threshold)

gap = diff(featureSet(:,1));

for i = numel(gap):-1:1
    if i <= indexUpperLimit
        break
    end
end

indexWidestGap = i; %index widest gap south of the algo determined threshold

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
    plot ([widestGapThreshold widestGapThreshold], ylim); %plot vertical line at 1
end

end

