function [events, thresholdFrequency, thresholdIntensity, thresholdDuration, indexArtifact, thresholdAmplitudeOutlier] = classifier_in_vitro (events, plotFigures)
%classifier_In_Vitro is hardcoded to detect ictal events
%   Classifying SLEs with hard-coded thresholds, very useful if no IIEs are
%   present to use the generic classifier.

%% Set Variables according to Michael's terms
userInput(3) = plotFigures;  %1 = yes; 0 = no
indexEventsToAnalyze = events(:,7)<4;

if nargin < 2
    userInput(3) = 0    %by default don't plot any figures
end

%% Remove artifacts (outliers) based on Peak-to-Peak Amplitude 
featureSet = events(:,6);   %peak-to-peak amplitude values
thresholdAmplitudeOutlier = 4.9;  %Michael's custom threshold
indexArtifact = events(:,6) > thresholdAmplitudeOutlier;  
index = indexArtifact; %Generic Terms

if sum(indexArtifact)>0   
    %Plot figure if artifacts detected within events
    if userInput(3) == 1      
        figArtifact = figure;
        gscatter(events(:,6) , events(:,6), index);    %plot index determined by Michael's Threshold
        hold on
        %plot Michael Chang's threshold values 
        plot ([thresholdAmplitudeOutlier thresholdAmplitudeOutlier], ylim);
        %Label
        set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
        set(gcf,'Name', 'Remove events containing artifact, using Peak-to-Peak Amplitude (mV)'); %select the name you want
        title ('Unsupervised classication, using k-means clustering');
        ylabel ('Peak-to-Peak Amplitude (mV)');
        xlabel ('Peak-to-Peak Amplitude (mV)');   
        legend('Epileptiform Events', 'Artifact', 'Michaels Artifact Threshold')
        legend ('Location', 'southeast')
        set(gca,'fontsize',12)        
    end

    %Remove artifact, based on Michael's threshold 
    events (indexArtifact, 7) = 4;  %%Label the event containing an artifact as '4'
    events (indexArtifact, 12) = 1;
    %Make new index without the artifacts
    indexEventsToAnalyze = events(:,7)<4;
end
    
%% classify based on frequency 
featureSet = events(:,4);   %Average Spike Rate (Hz)
%Michael's threshold
michaelsFrequencyThreshold = 1; %Hz  
%Algo determined threshold
[algoFrequencyIndex, algoFrequencyThreshold] = findThresholdSLE (events(indexEventsToAnalyze,4));
%Use Michael's hard-coded threshold 
thresholdFrequency = michaelsFrequencyThreshold;    %michael's threshold frequency is the lowest frequency for SLEs

%Event is a SLE if larger than threshold
indexFrequency = featureSet>=thresholdFrequency;  
events (:,9) = indexFrequency;
%Plot figure
    if userInput(3) == 1  
    %plot figure
    index = indexFrequency;   %convert to generic name for plotting
    featureThreshold = thresholdFrequency;  %convert to generic name for plotting
    figFrequency = figure;
    gscatter(featureSet , featureSet, index);      
    hold on
    %plot the algorithm detected threshold
    plot ([algoFrequencyThreshold algoFrequencyThreshold], ylim); 
    %plot Michael Chang's threshold  
    plot ([michaelsFrequencyThreshold michaelsFrequencyThreshold], ylim);
    %plot threshold that was used
    plot ([featureThreshold featureThreshold], ylim);

    %Label
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Feature Set: Spiking Rate (Hz)'); %select the name you want    
    title ('Unsupervised classication, using k-means clustering');
    ylabel ('Spiking Rate (Hz)');
    xlabel ('Spiking Rate (Hz)');  
        if numel(unique(indexFrequency))>1  %legend depends on what's present
            legend('IIE', 'SLE', 'Algo Threshold', 'Michaels Threshold', 'Frequency Threshold')
        else if unique(indexFrequency)==1
                legend('SLE', 'Algo Threshold', 'Michaels Threshold', 'Frequency Threshold')
            else
                legend('IIE', 'Algo Threshold', 'Michaels Threshold', 'Frequency Threshold')
            end
        end                                  
    set(gca,'fontsize',12)
    end

%% classify based on intensity 
featureSet = events(:,5);   %Average intensity (Power/Duration)
%Michael's threshold
if mean(events(indexEventsToAnalyze,5))>std(events(indexEventsToAnalyze,5))
    michaelIntensityThreshold = mean(events(indexEventsToAnalyze,5))-std(events(indexEventsToAnalyze,5));
else 
    michaelIntensityThreshold = mean(events(indexEventsToAnalyze,5));
end
%Algo determined threshold 
[algoIntensityIndex, algoIntensityThreshold] = findThresholdSLE (events(indexEventsToAnalyze,5));

%use the hard-coded threshold for Intensity, (floor: 5 mV^2/s)
thresholdIntensity = 5;

%determine the index for SLE and IIE using threshold for Intensity (feature)
indexIntensity = featureSet>=thresholdIntensity;
events (:,10) = indexIntensity; %store in array

    if userInput(3) == 1  
    %plot figure
    index = indexIntensity;
    featureThreshold = thresholdIntensity;
    figIntensity = figure;
    gscatter(featureSet , featureSet, index);    %plot scatter plot
    hold on
    %plot the algorithm detected threshold
    plot ([algoIntensityThreshold algoIntensityThreshold], ylim); 
    %plot Michael Chang's threshold values 
    plot ([michaelIntensityThreshold michaelIntensityThreshold], ylim);
    %plot threshold that was used
    plot ([featureThreshold featureThreshold], ylim);

    %Label
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Feature Set: Intensity (Power/Duration)'); %select the name you want   
    title ('Unsupervised classication, using k-means clustering');
    ylabel ('Average Intensity (Power/Duration)');
    xlabel ('Average Intensity (Power/Duration)');   
        if numel(unique(indexIntensity))>1  %legend depends on what's present
            legend('IIE', 'SLE', 'Algo Threshold', 'Michaels Threshold', 'Intensity Threshold')
        else if unique(indexIntensity)==1
                legend('SLE', 'Algo Threshold', 'Michaels Threshold', 'Intensity Threshold')
            else
                legend('IIE', 'Algo Threshold', 'Michaels Threshold', 'Intensity Threshold')
            end
        end    
    set(gca,'fontsize',12)
    end

%% Classify based on duration 
featureSet = events(:,3);   %Duration (s)

%Average epileptiform event duration and sigma
averageSLEDuration = mean(featureSet);
sigmaSLEDuration = std(featureSet);

%Michael's threshold, use the one that is higher, conservative
if averageSLEDuration-(2*sigmaSLEDuration) > sigmaSLEDuration
    michaelsDurationThreshold=averageSLEDuration-(2*sigmaSLEDuration);
else
    michaelsDurationThreshold=sigmaSLEDuration;  
end

%Algo deteremined threshold (tend to be higher value)
[algoDurationIndex, algoDurationThreshold] = findThresholdSLE (events(indexEventsToAnalyze,3));

%Use the hard-coded threhsold, 10 s (the floor)
thresholdDuration = 10;

indexDuration = featureSet>thresholdDuration; 
events(:,11) = indexDuration;

    if userInput(3) == 1  
    %plot figure
    index = indexDuration;
    featureThreshold = thresholdDuration;
    figDuration = figure;
    gscatter(featureSet , featureSet, index);    %plot scatter plot
    hold on
    %plot the algorithm detected threshold
    plot ([algoDurationThreshold algoDurationThreshold], ylim); 
    %plot Michael Chang's threshold values 
    plot ([michaelsDurationThreshold michaelsDurationThreshold], ylim);
    %plot threshold that was used
    plot ([featureThreshold featureThreshold], ylim); 
    %Label
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', 'Feature set: Duration (sec)'); %select the name you want
    title ('Unsupervised classication, using k-means clustering');
    ylabel ('Duration (sec)');
    xlabel ('Duration (sec)');   
        if numel(unique(indexDuration))>1  %legend depends on what's present
            legend('IIE', 'SLE', 'Algo Threshold', 'Michaels Threshold', 'Duration Threshold')
        else if unique(indexDuration)==1
                legend('SLE', 'Algo Threshold', 'Michaels Threshold', 'Duration Threshold')
            else
                legend('IIE', 'Algo Threshold', 'Michaels Threshold', 'Duration Threshold')
            end
        end 
    set(gca,'fontsize',12)
    end

%% Classification (final)
for i = 1: numel(events(:,1))
    if indexFrequency(i) + indexIntensity(i) + indexDuration(i) == 3 & events (i,7) <4
        events (i,7) = 1;   %1 = SLE; 2 = IIE; 3 = IIS; 0 = unclassified.
        else if indexFrequency(i) + indexIntensity(i) + indexDuration(i) < 3 & events (i,7) < 4
                events (i,7) = 2;
                else
                events (i,7) = 4;
            end
    end
end


end

