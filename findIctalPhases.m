function [eventPhases,spikeFrequency] = findIctalPhases(spikeFrequency,frequency)
%findIctalPhases identifies the different phases of an epileptiform event.
%   The ictal event is comprised of the sentinel spike, the preictal phase, 
%   tonic-like phase, and clonic-like phase. The tonic-like firing is the 
%   period of the ictal event where there is contiguous frequency (>2 s)
%   above or equal to 1/3 maximum frequency of the event. The tonic phase
%   is considered over when there is continguous frequency (>2 s) below 1/3
%   maximum frequency of the event or if frequency is 0 for at least 1 s.
%   The clonic phase is characterized as the period following the tonic
%   phase when the frequency is inconsistent interspersed with values of 0
%   throughout the phase. The preictal phase is the period between the
%   sentinel spike and start of the tonic-phase. This function uses the
%   frequency feature set (spike rate per second) to chracterize the
%   different phases of an ictal event. If no tonic phase is identified 
%   The intensity feature set can also reveal additional insights into the 
%   shape and morphology of the ictal event. However, it is currently not
%   implemented. Additional Notes about the output:
%     eventPhases(1) = startTonicTime;
%     eventPhases(2) = endTonicTime;
%     eventPhases(3) = classification;
%     eventPhases(4) = preictalPhaseDuration;
%     eventPhases(5) = tonicPhaseDuration;
%     eventPhases(6) = clonicPhaseDuration;
%     eventPhases(7) = startTime;
%     eventPhases(8) = endTime;

%% Set default values if not specified
if nargin <2
    frequency = 10000;  %Hz
end

if numel(spikeFrequency(:,1)) > 3
    %% Use freqency feature set to classify epileptiform event
    maxFrequency = double(max(spikeFrequency(:,2)));    %Calculate Maximum frequency during the SLE
    indexTonic = spikeFrequency(:,2) >= maxFrequency/3; %Use Michael's threshold to seperate frequency feature set into two populations, high and low.
    spikeFrequency(:,3) = indexTonic;    %store Boolean index 

    %locate start of Tonic phase | Contingous segments above threshold    
    for j = 2: numel (indexTonic) %slide along the SLE; ignore the first spike, which is the sentinel spike
        if j+1 > numel(indexTonic)  %If you scan through the entire SLE and don't find a tonic phase, reclassify event as a IIE            
            j = find(indexTonic(2:end),1,'first');  %Take the first second frequency is 'high' as onset if back-to-back high frequency are not found
                if isempty(j)
                    j = 1;  %honory position, just to push the function through
                end                
            startTonic = spikeFrequency(j+1);  %j+1 because you started the find function (line above) starting from position "2:end"
            endTonic = spikeFrequency(numel(indexTonic));    %if no tonic period is found; just state whole ictal event as a tonic phase
            classification = 0; % 0 = no tonic phase, 1 = tonic-clonic SLE, 2 - tonic-only        
        else                        
            if indexTonic(j) > 0 && indexTonic(j+1) > 0 %If you locate two continuous segments with high frequency, mark the first segment as start of tonic phase                        
                startTonic = spikeFrequency(j);  %store the onset time
                classification = 1;   %1 = tonic-clonic SLE;   2 = tonic-only
                while ~and(indexTonic(j) == 0, indexTonic(j+1) == 0) & spikeFrequency(j,2) ~= 0; %Locate the offset time as the first point as either two continueous segments with low frequency or zero frequency, which ever comes first
                    j = j+1;    %keep sliding along the SLE until the statement above is false.
                    if j+1 > numel(indexTonic)  %If you slide all the way to the end and still can't find tonic phase offset, 
                        j = numel(indexTonic)+1;  %take the last point as the offset of the tonic phase - this means there is no clonic phase; add 1 because note you will remove it in line 33
                        classification = 2;   %1 = SLE;   2 = Tonic-only     
                        break
                    end                                    
                end            
                endTonic = spikeFrequency(j-1);
                break
            end
        end        
    end          

    %Characterize the epileptiform event | Label all the features of the SLE
    startTime = spikeFrequency(1,1)/frequency;  %secs
    endTime = spikeFrequency(end,1)/frequency;  %secs
    startTonicTime = startTonic/frequency;  %secs
    endTonicTime = endTonic/frequency;  %secs

    preictalPhaseDuration = startTonicTime - startTime;
    tonicPhaseDuration = endTonicTime  - startTonicTime;
    clonicPhaseDuration = endTime - endTonicTime; 
else
    startTime = spikeFrequency(1,1)/frequency;  %secs
    endTime = spikeFrequency(end,1)/frequency;  %secs
    startTonicTime = startTime;  %secs
    endTonicTime = endTime;  %secs
    
    classification = -1;    %-1 = IIS (<3s duration)  
    fprintf(2,'\nWarning: The event at %d-%d s is not a IIE (duration <2 s) so it was reclassified as a IIS (-1).\n', startTime, endTime)

    preictalPhaseDuration = startTonicTime - startTime;
    tonicPhaseDuration = endTonicTime  - startTonicTime;
    clonicPhaseDuration = endTime - endTonicTime; 
end
    
%Store event's characteristics for output
eventPhases(1) = startTonicTime;
eventPhases(2) = endTonicTime;
eventPhases(3) = classification;
eventPhases(4) = preictalPhaseDuration;
eventPhases(5) = tonicPhaseDuration;
eventPhases(6) = clonicPhaseDuration;
% eventPhases(7) = startTime;
% eventPhases(8) = endTime;



