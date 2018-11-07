function [thresholdAmplitudeOutlier ] = findThresholdArtifact(featureSet,plotGraph)
%[indexEvents,widestGapThreshold ] = findThresholdArtifact(featureSet,plotGraph)
%   This function classifies a populations of events into two distinct
%   populations based one their peak-to-peak amplitude feature set. This 
%   function detects if any outliers exist 3.2 standard deviations above 
%   the average value. If an outlier is detected, it will find the
%   threshold as the mid-point of the widest gap below the lowest value
%   outlier. This method is to ensure outliers close to the artifact 
%   threshold are not missed and accounted for because there is always a
%   slight margin of error (give) that exists. The function will sort the 
%   feature set and seperate the population based on where the widest gap
%   (in the value of the feature set) is observed, at or below the outlier
%   with the lowest value. The population of events with a larger
%   value (of the feature) is classied an artifact (outlier), and the
%   population with a lower value is classied as epileptiform events. 
%   There is a feature to plot %   a graph of the feature set and the 
%   thresholds found by this function algorithm. Set the second input 
%   variable to '1' to plot graph for troubleshooting.


if nargin == 1  
    plotGraph = 0;  %1 = yes; 0 = no
end   

%Michael's custom threshold, to detect if outlier is present
michaelAmplitudeThreshold = mean(featureSet)+(3.2*std(featureSet));  %Michael's custom threshold
michaelAmplitudeIndex = featureSet > michaelAmplitudeThreshold;  
thresholdAmplitudeOutlier = michaelAmplitudeThreshold;

if sum(michaelAmplitudeIndex)>0    %check if outlier is present
    %% Find widest gap (below the detected outlier)
    featureSet = sort(featureSet);   
   
    %find index of upper limit using Michael's threshold
    indexUpperLimit = find(featureSet>michaelAmplitudeThreshold, 1); 
    
    gap = diff(featureSet);    %find the size of gaps
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
    thresholdAmplitudeOutlier = widestGapThreshold; 
end

%segment events into Artifacts and Epileptiform Events, using widest gap threshold
indexEvents = featureSet>thresholdAmplitudeOutlier;

%% plot figure

if plotGraph == 1    
    figure;
    gscatter(featureSet , featureSet, indexEvents);    %plot scatter plot
    hold on
    plot ([michaelAmplitudeThreshold michaelAmplitudeThreshold], ylim);    
    plot ([thresholdAmplitudeOutlier thresholdAmplitudeOutlier], ylim); 
    legend ('Artifacts', 'Epileptiform Event', 'Michaels Threshold', 'widest-Gap Threshold')
    title ('Unsupervised Artifact threshold for classification')
end

end

