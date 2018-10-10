function [events, thresholdFrequency, thresholdIntensity, thresholdDuration, indexArtifact, thresholdAmplitudeOutlier, algoFrequencyThreshold] = classifier_dynamic (events, plotFigures, IIE_classifier, floorThresholdFrequency, hardCodedThreshold)
%classifier_dynamic will classify a population (>6) of epilpeptiform events
%into two different populations (IIEs and SLEs)
%   The default setting for this function is to classify ictal events as
%   SLEs. The floor threshold for frequency is 1 Hz. If you switch it to
%   classify IIEs, the floor threshold for potential SLEs will drop to 0.6
%   Hz.
%   The threshold to classify populations of epileptiform events are based
%   are drawn from three different methods of calculation, hard-codes based
%   on the in vitro ictal event population, Michael's calculation based on
%   emperical data, and k-means clustering (in combinations with the widest
%   gap). Note: intensity and duration do not require any floors, they tend
%   to balance each other out, however, the floor threshold for intensity
%   is set to 1 mV^2/s. I put this in place incase it was analyzing only
%   interictal spikes and events.

%% Set Variables according to Michael's terms
userInput(3) = plotFigures;  %1 = yes; 0 = no

%% Set default variables, if not specified
if nargin < 5    
    hardCodedThreshold = numel(events(:,1));    %turned on if less than 1; >=2 it is Off
end

if nargin < 4    
    IIE_classifier = 0; %Turned off (classify SLEs)
    hardCodedThreshold = numel(events(:,1));   
end

if nargin < 3
    userInput(3) = 0;    %by default don't plot any figures
    IIE_classifier = 0; 
    hardCodedThreshold = numel(events(:,1));
end

%% Switch between classifier for SLEs and IIEs
switch IIE_classifier            
    case 1
        indexEventsToAnalyze = events(:,7)==2;  %Analyze IIEs
        %Indicate any issues
        if hardCodedThreshold < 1
            fprintf(2,'\nNo interictal events were detected for further analysis.\n')
            return
        elseif hardCodedThreshold < 2
            fprintf(2,'\nOnly 1 interictal event detected; interictal event(s) will be classified using hard-coded thresholds based on Chang et al., 2018 Neurobiology of Disease.\n')                    
        end
        %Set Floor threshold, if not specified
        if nargin < 4
            floorThresholdFrequency = 0.6;    %by default don't plot any figures            
        end
        if nargin < 3            
            floorThresholdFrequency = 0.6;    %Hz            
        end
    
    case 0
        indexEventsToAnalyze = events(:,7)<4;   %Analyze SLEs
        %Indicate any issues      
        if hardCodedThreshold < 1
            fprintf(2,'\nNo epileptiform events were detected.\n')
            return
        elseif hardCodedThreshold < 2
            fprintf(2,'\nOnly 1 events was detected in the in vitro recording; the events will be classified using hard-coded thresholds based on Chang et al., 2018 Neurobiology of Disease.\n')        
        end
        %Set Floor threshold, if not specified
        if nargin < 4
            floorThresholdFrequency = 1;    %by default don't plot any figures            
        end
        if nargin < 3            
            floorThresholdFrequency = 1;    %Hz            
        end

    %% Stage 1: Artifact (Amplitude Outlier) removal | only for SLE Classifier
    %Remove outliers based on peak-to-peak amplitude
    featureSet = events(:,6);   %peak-to-peak amplitude values

    %algo-determined threshold | Find widest gap (below the detected outlier)
    if hardCodedThreshold>=2
        thresholdAmplitudeOutlier = findThresholdArtifact (featureSet);
    else
        thresholdAmplitudeOutlier = 4.9;    %Michael's custom threshold
    end
    
    %segment events into Artifacts and Epileptiform Events, using algo-deteremined threshold
    indexArtifact = featureSet>thresholdAmplitudeOutlier;
    index = indexArtifact; %Generic Terms

    if sum(index)>0   
        %Plot figure if artifacts detected within events
        if userInput(3) == 1      
            figArtifact = figure;
            gscatter(events(:,6) , events(:,6), index);    %plot index determined by Michael's Threshold
            hold on
            %plot the determined threshold values 
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

        %Remove artifact
        events (indexArtifact, 7) = 4;  %%Label the event containing an artifact as '4'
        events (indexArtifact, 12) = 1;
        %Make new index without the artifacts
        indexEventsToAnalyze = events(:,7)<4;
    end
end

%% Stage 2: Unsupervised Classifier 
% classify based on average frequency 
featureSet = events(:,4);   %Average Spike Rate (Hz)

%Find threshold, dynamic 
michaelsFrequencyThreshold = floorThresholdFrequency;   %Michael's threshold (floor threshold)
if hardCodedThreshold>= 2
    [~, algoFrequencyThreshold] = findThresholdSLE (events(indexEventsToAnalyze,4));    %Algo determined threshold; use this as the floor when spliting IIEs
else
    algoFrequencyThreshold = floorThresholdFrequency;    %If it is unable to perform the calculate, just use the floor threshold.
end

%Use the largest threshold
thresholdFrequency = max([michaelsFrequencyThreshold, algoFrequencyThreshold]);
    
%Event is a SLE if larger than threshold
indexFrequency = featureSet>=thresholdFrequency;    %split the population based on feature  
events (:,9) = indexFrequency;  %store the index

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
    legend ('Location', 'southeast')
    set(gca,'fontsize',12)
    end
    
    
%classify based on average intensity
featureSet = events(:,5);   %Average intensity (Power/Duration)
%Floor threshold
floorIntensityThreshold = 5;    %mV^2/s
%Dynamic threshold
if hardCodedThreshold >=2
    %Michael's threshold
    if mean(events(indexEventsToAnalyze,5))>std(events(indexEventsToAnalyze,5))
        michaelIntensityThreshold = mean(events(indexEventsToAnalyze,5))-std(events(indexEventsToAnalyze,5));
    else 
        michaelIntensityThreshold = mean(events(indexEventsToAnalyze,5));
    end
    %Algo determined threshold
    [~, algoIntensityThreshold] = findThresholdSLE (events(indexEventsToAnalyze,5));
else
    michaelIntensityThreshold = floorIntensityThreshold;
    algoIntensityThreshold = michaelIntensityThreshold;
end
    
%use the lower threshold for Intensity, (with a floor threshold in place at 5 mV^2/s)
if algoIntensityThreshold <= floorIntensityThreshold || michaelIntensityThreshold <= floorIntensityThreshold
    thresholdIntensity = floorIntensityThreshold;
else if algoIntensityThreshold<=michaelIntensityThreshold
        thresholdIntensity = algoIntensityThreshold;
    else
        thresholdIntensity = michaelIntensityThreshold;
    end
end

%determine the index for SLE and IIE using threshold for Intensity (feature)
indexIntensity = featureSet>=thresholdIntensity;    %split the population based on the feature
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
          
    legend ('Location', 'southeast')
    set(gca,'fontsize',12)
    end

%% Switch between SLE classifier and IIE classifier
switch IIE_classifier 
    case 0
        %Initial SLE Classifier (Filter)
        for i = 1: numel(events(:,1))
            if (indexFrequency(i) + indexIntensity(i)) == 2 & events(i,7) <4    %and not an artifact
                events (i,7) = 1;   %1 = SLE; 2 = IIE; 3 = IIS; 0 = unclassified
            else if (indexFrequency(i) + indexIntensity(i)) < 2 & events(i,7) <4
                    events (i,7) = 2;   %2 = IIE
                else            
                    events (i,7) = 4;   %4 = event contains an artifact
                end
            end    
        end

    case 1
        %Initial IIE Classifier (Filter)
        for i = 1: numel(events(:,1))
            if (indexFrequency(i) + indexIntensity(i)) == 2 & events(i,7) <4    %and not an artifact
                events (i,7) = 1;   %1 = SLE; 2 = IIE; 3 = IIS; 0 = unclassified; 4 = artifact; 5 = Class 5 SLE (questionable)
            else
                events (i,7) = 2;   %2 = IIE
            end        
        end  
end


%% Stage 3: Machine Learning to classify using duration as a feature
%Putative SLEs
indexPutativeSLE = events (:,7) == 1;  %classify which ones are SLEs
putativeSLE = events(indexPutativeSLE, :);

%Find threshold for duration, using Michael's method (form of Machine Learning)
averageSLEDuration=mean(putativeSLE(:,3));
sigmaSLEDuration = std(putativeSLE(:,3));

%classify based on duration
featureSet = events(:,3);   %Duration (s)
%Floor threshold
floorThresholdDuration = 10;    %seconds
%Dynamic threshold
if hardCodedThreshold>=2
    %Michael's threshold, use the one that is higher, conservative
    if averageSLEDuration-(2*sigmaSLEDuration) > sigmaSLEDuration
        michaelsDurationThreshold=averageSLEDuration-(2*sigmaSLEDuration);
    else
        michaelsDurationThreshold=sigmaSLEDuration;  
    end
    %Algo deteremined threshold (tend to be higher value)
    [~, algoDurationThreshold] = findThresholdSLE (events(indexEventsToAnalyze,3));
else
    michaelsDurationThreshold = floorThresholdDuration;
    algoDurationThreshold = michaelsDurationThreshold;
end

%Use the lowest (more liberal) threhsold, unless it's below (the floor)
if algoDurationThreshold <= floorThresholdDuration || michaelsDurationThreshold <= floorThresholdDuration
    thresholdDuration = floorThresholdDuration;
else if michaelsDurationThreshold <= algoDurationThreshold
        thresholdDuration = michaelsDurationThreshold;
    else
        thresholdDuration = algoDurationThreshold;
    end
end

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
    legend ('Location', 'southeast')
    set(gca,'fontsize',12)
    end

%% Final Classifier | Swtich between SLE classifier and IIE classifier
switch IIE_classifier
    case 0
        %Final SLE Classifier
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

    case 1
        %Final IIE Classifier
        for i = 1: numel(events(:,1))
            if indexFrequency(i) + indexIntensity(i) + indexDuration(i) == 3 & events (i,7) <4
                events (i,7) = 1.5;   %1.5 = Questionable seizure; 2.5 = Questionable IIE; 0 = unclassified; 
            else
                events (i,7) = 2.5;
            end
        end
end 



