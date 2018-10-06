%% Confirm class 1 SLEs, by locating tonic-like and clonic-like firing
% for i = indexSLE
%     [tonicPhase, spikeFrequency{i}] = findIctalPhases (spikeFrequency{i});
% end
             
  %% Locate SLE phases
    for i = 1:numel(events(:,1))
  %Use freqency feature set to classify SLEs
    maxFrequency = double(max(spikeFrequency{i}(:,2)));    %Calculate Maximum frequency during the SLE
    indexTonic = spikeFrequency{i}(:,2) >= maxFrequency/3; %Use Michael's threshold to seperate frequency feature set into two populations, high and low.
    spikeFrequency{i}(:,3) = indexTonic;    %store Boolean index 

    %locate start of Tonic phase | Contingous segments above threshold    
    for j = 2: numel (indexTonic) %slide along the SLE; ignore the first spike, which is the sentinel spike
        if j+1 > numel(indexTonic)  %If you scan through the entire SLE and don't find a tonic phase, reclassify event as a IIE            
            j = find(indexTonic(2:end),1,'first');  %Take the first second frequency is 'high' as onset if back-to-back high frequency are not found
            startTonic = spikeFrequency{i}(j);  %store the onset time                                   
            endTonic = spikeFrequency{i}(numel(indexTonic));    %if no tonic period is found; just state whole ictal event as a tonic phase
            events (i,7) = 2;   %Update the event's classification as a IIE (since its lacking a tonic phase)
            events (i,13) = 0;
            fprintf(2,'\nWarning: The Tonic Phase was not detected in SLE #%d from File: %s; so it was reclassified as a IIE.\n', i, FileName)
        else                        
            if indexTonic(j) > 0 && indexTonic(j+1) > 0 %If you locate two continuous segments with high frequency, mark the first segment as start of tonic phase                        
                startTonic = spikeFrequency{i}(j);  %store the onset time
                events(i,13) = 1;   %1 = tonic-clonic SLE;   2 = tonic-only
                while ~and(indexTonic(j) == 0, indexTonic(j+1) == 0) & spikeFrequency{i}(j,2) ~= 0; %Locate the offset time as the first point as either two continueous segments with low frequency or zero frequency, which ever comes first
                    j = j+1;    %keep sliding along the SLE until the statement above is false.
                    if j+1 > numel(indexTonic)  %If you slide all the way to the end and still can't find tonic phase offset, 
                        j = numel(indexTonic)+1;  %take the last point as the offset of the tonic phase - this means there is no clonic phase
%                         endTonic = spikeFrequency{i}(j);    %if no tonic period is found; just state whole ictal event as a tonic phase
                        events(i,13) = 2;   %1 = SLE;   2 = Tonic-only                        
                        fprintf(2,'\nWarning: No clonic period detected in SLE #%d from File: %s. Review your data.\n', i, FileName)                        
                        break
                    end                                    
                end            
                endTonic = spikeFrequency{i}(j-1);
                break
            end
        end        
    end        

%     %Label all the features of the SLE
%     startSLETime = spikeFrequency{i}(1,1);
%     endSLETime = spikeFrequency{i}(end,1);
%     startTonicTime = startTonic/frequency;
%     endTonicTime = endTonic/frequency;
%     
%     preictalPhase = startTonicTime - startSLETime;
%     tonicPhase = endTonicTime  - startTonicTime;
%     clonicPhase = endSLETime - endTonicTime; 
%     
%     events(i,13) = spikeFrequency{i}(end,1)
%     events(i,14) =
%     events(i,15) = 
    
    if userInput(5) == 1     
    %% Using k-means clustering (algo threshold) for intensity | Store for future use (when you know what to do with it)
%     k = 2;
%     maxIntensity = double(max(intensityPerMinute{i}(:,2)));
%     indexAnalyze = intensityPerMinute{i}(:,2) > (maxIntensity/10); 
%     featureSet = intensityPerMinute{i}(indexAnalyze,2);
%     indexIntensityAlgo = kmeans(featureSet, k);
%     intensityPerMinute{i}(:,4) = 0; 
%     intensityPerMinute{i}(indexAnalyze,4) = indexIntensityAlgo; 
    
% %     %locate contingous segments above threshold    
%     for i = 2: numel (indexTonic)
%         if indexTonic(i) > 0 && indexTonic(i+1) > 0                        
%             startTonic = spikeFrequency{i}(i);
%             while indexTonic(i) > 0
%                 i = i+1;
%             end            
%             endTonic = spikeFrequency{i}(i-1);
%             break
%         end
%     end       
    
    %Plot Figure
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('Frequency Feature Set, Epileptiform Event #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize')); 

%     subplot (2,1,1)
    figHandle = plotEvent (figHandle, LFP_centered, t, events(i,1:2), locs_spike_2nd, lightpulse);        
    %Labels
    title (sprintf('Epileptiform Event #%d, Michaels Threshold', i));
    ylabel ('mV');
    xlabel ('Time (sec)');   
    %Plot Frequency Feature
    yyaxis right
    
    activeIndex = spikeFrequency{i}(:,3) == 1;
    inactiveIndex = spikeFrequency{i}(:,3) == 0;
    
    plot (spikeFrequency{i}(:,1)/frequency, spikeFrequency{i}(:,2), 'o', 'MarkerFaceColor', 'cyan')
    
    plot (spikeFrequency{i}(inactiveIndex ,1)/frequency, spikeFrequency{i}(inactiveIndex ,2), 'o', 'MarkerFaceColor', 'magenta')
    plot (spikeFrequency{i}(inactiveIndex ,1)/frequency, spikeFrequency{i}(inactiveIndex ,2), 'o', 'color', 'k')
    
    plot ([startTonic/frequency startTonic/frequency], ylim)
    plot ([endTonic/frequency endTonic/frequency], ylim)
    
    ylabel ('Spike rate/second (Hz)');
    set(gca,'fontsize',14)
    legend ('LFP filtered', 'Epileptiform Event', 'Detected Onset', 'Detected Offset', 'Detected Spikes', 'Applied Stimulus', 'Frequency above threshold', 'Frequency below threshold')
    legend ('Location', 'northeastoutside')

    
%     subplot (2,1,2)    
%     figHandle = plotEvent (figHandle, LFP_centered, t, events(i,1:2), locs_spike_2nd, lightpulse);        
%     %Labels
%     title ('Algo Threshold (K-means clustering)');
%     ylabel ('mV');
%     xlabel ('Time (sec)');   
%     %Plot Frequency Feature
%     yyaxis right
%     
%     fourthIndex = intensityPerMinute{i}(:,4) == 4;
%     bottomIndex = intensityPerMinute{i}(:,4) == 3;
%     middleIndex = intensityPerMinute{i}(:,4) == 2;
%     topIndex = intensityPerMinute{i}(:,4) == 1;
%     zeroIndex = intensityPerMinute{i}(:,4) == 0;
% 
%     plot (intensityPerMinute{i}(:,1)/frequency, intensityPerMinute{i}(:,2), 'o', 'MarkerFaceColor', 'cyan')    
%     plot (intensityPerMinute{i}(middleIndex ,1)/frequency, intensityPerMinute{i}(middleIndex ,2), 'o', 'MarkerFaceColor', 'yellow')    
%     plot (intensityPerMinute{i}(bottomIndex ,1)/frequency, intensityPerMinute{i}(bottomIndex ,2), 'o', 'MarkerFaceColor', 'magenta')
%     plot (intensityPerMinute{i}(fourthIndex ,1)/frequency, intensityPerMinute{i}(fourthIndex ,2), 'o', 'MarkerFaceColor', 'black')
%     plot (intensityPerMinute{i}(zeroIndex ,1)/frequency, intensityPerMinute{i}(zeroIndex ,2), 'o', 'MarkerFaceColor', 'black')
% 
%     plot (intensityPerMinute{i}(:,1)/frequency, intensityPerMinute{i}(:,2), 'o', 'color', 'k')
%     
%     ylabel ('Intensity (mW^2/s)');
%     set(gca,'fontsize',14)
%     legend ('LFP filtered', 'Epileptiform Event', 'Detected Onset', 'Detected Offset', 'Detected Spikes', 'Applied Stimulus', 'Intensity - Gp1', 'Intensity - Gp2', 'Intensity - Gp3')
%     legend ('Location', 'northeastoutside')

    %Export figures to .pptx
    exportToPPTX('addslide'); %Draw seizure figure on new powerpoint slide
    exportToPPTX('addpicture',figHandle);      
    close(figHandle)
    end
    end
       
if userInput(5) == 1
    % save and close the .PPTX
    exportToPPTX('saveandclose',sprintf('%s(freq,intensity)', excelFileName)); 
end