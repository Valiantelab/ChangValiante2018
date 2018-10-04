function [tonicPhase, spikeRateMinute] = findIctalPhases(spikeRateMinute,frequency)
%findIctalPhases function identifies the different ictal event phases
%   The ictal event is comprised of the sentinel spike, the preictal phase, 
%   tonic-like phase, and clonic-like phase. The tonic-like firing is the 
%   period of the ictal event where there is contiguous frequency 
%   above the threshold of 1/3 maximum frequency of the ictal event. The
%   tonic-like phase is characterized as the period following the tonic
%   phase when the frequency is inconsistent intersperced with values of 0
%   throughout the phase. The preictal phase is the period of time between the
%   sentinel spike and start of the tonic-phase. This function uses the
%   frequency feature set (spike rate per second) to chracterize the
%   different phases of an ictal event. The intensity feature set can
%   reveal additional insights into the shape and morphology of the ictal
%   event.

%% Default values if not specified
if nargin <2
    frequency = 10000;
end

        
    %% Using Michael's Threshold    
    %locate putatitve tonic phase
    maxFrequency = double(max(spikeRateMinute(:,2)));    %Calculate the maximum frequency         
    indexTonic = spikeRateMinute(:,2) > maxFrequency/3; %Locate putative Tonic Phase    
    spikeRateMinute(:,3) = indexTonic;
    
    %locate contingous segments above threshold    
    for i = 2: numel (indexTonic)
        if indexTonic(i) > 0 && indexTonic(i+1) > 0                        
            startTonicTime = spikeRateMinute(i);
            while indexTonic(i) > 0
                i = i+1;
            end            
            endTonicTime = spikeRateMinute(i-1);
            break
        end
    end          
    
    tonicPhase (1,1) = startTonicTime;
    tonicPhase (1,2) = endTonicTime;
    
    %% Using k-means clustering (algo threshold)
    k = 3;
    featureSet = spikeRateMinute(:,2);
    indexFrequencyAlgo = kmeans(featureSet, k);
    spikeRateMinute(:,4) = indexFrequencyAlgo; 
    
    %locate contingous segments above threshold    
    for i = 2: numel (indexTonic)
        if indexTonic(i) > 0 && indexTonic(i+1) > 0                        
            startTonicTime = spikeRateMinute(i);
            while indexTonic(i) > 0
                i = i+1;
            end            
            endTonicTime = spikeRateMinute(i-1);
            break
        end
    end       
    
    %Plot Figure
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('Frequency Feature Set, Epileptiform Event #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize')); 

    subplot (2,1,1)
    figHandle = plotEvent (figHandle, LFP_centered, t, events(i,1:2), locs_spike_2nd, lightpulse);        
    %Labels
    title (sprintf('Epileptiform Event #%d, Michaels Threshold', i));
    ylabel ('mV');
    xlabel ('Time (sec)');   
    %Plot Frequency Feature
    yyaxis right
    
    activeIndex = spikeRateMinute(:,3) == 1;
    inactiveIndex = spikeRateMinute(:,3) == 0;
    
    plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'MarkerFaceColor', 'cyan')
    
    plot (spikeRateMinute(inactiveIndex ,1)/frequency, spikeRateMinute(inactiveIndex ,2), 'o', 'MarkerFaceColor', 'magenta')
    plot (spikeRateMinute(inactiveIndex ,1)/frequency, spikeRateMinute(inactiveIndex ,2), 'o', 'color', 'k')
    
    plot ([startTonicTime/frequency startTonicTime/frequency], ylim)
    plot ([endTonicTime/frequency endTonicTime/frequency], ylim)
    
    ylabel ('Spike rate/second (Hz)');
    set(gca,'fontsize',14)
    legend ('LFP filtered', 'Epileptiform Event', 'Detected Onset', 'Detected Offset', 'Detected Spikes', 'Applied Stimulus', 'Frequency above threshold', 'Frequency below threshold')
    legend ('Location', 'northeastoutside')

    
    subplot (2,1,2)    
    figHandle = plotEvent (figHandle, LFP_centered, t, events(i,1:2), locs_spike_2nd, lightpulse);        
    %Labels
    title ('Algo Threshold (K-means clustering)');
    ylabel ('mV');
    xlabel ('Time (sec)');   
    %Plot Frequency Feature
    yyaxis right
    
    bottomIndex = spikeRateMinute(:,4) == 3;
    middleIndex = spikeRateMinute(:,4) == 2;
    topIndex = spikeRateMinute(:,4) == 1;
        
    plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'MarkerFaceColor', 'cyan')    
    plot (spikeRateMinute(middleIndex ,1)/frequency, spikeRateMinute(middleIndex ,2), 'o', 'MarkerFaceColor', 'yellow')    
    plot (spikeRateMinute(bottomIndex ,1)/frequency, spikeRateMinute(bottomIndex ,2), 'o', 'MarkerFaceColor', 'magenta')
    plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'color', 'k')
    
    ylabel ('Spike rate/second (Hz)');
    set(gca,'fontsize',14)
    legend ('LFP filtered', 'Epileptiform Event', 'Detected Onset', 'Detected Offset', 'Detected Spikes', 'Applied Stimulus', 'Frequency - Gp1', 'Frequency - Gp2', 'Frequency - Gp3')
    legend ('Location', 'northeastoutside')
 


%     %using intensity
%     maxIntensity = max(intensity)
%     
%     %Tonic-like firing
%     indexTonic = intensity(2:end,2) > maxIntensity/3
%     indexTonicStart = find(indexTonic(:,2),1,'first')
%     intensity(indexTonicStart+1,:)    %add 1 to add back in the sentinel spike (index) you removed
%     
%     %Clonic-like firing
%     indexClonic = intensity(2:end,2) > maxIntensity/1.33
%     indexClonicStart = find(indexClonic(:,2),1,'first')
%     intensity(indexTonicStart+1,:)    %add 1 to add back in the sentinel spike (index) you removed
%     
%     %k-means clustering
%     featureSet = intensity(:,2);
%      indexIntensity = kmeans(featureSet, 3);
%      
%      figure
%      gscatter(featureSet, featureSet, idx)

end

