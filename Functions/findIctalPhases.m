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

if numel(spikeFrequency(:,1)) > 3   %Analyze if the event is greater than 3 seconds
    %% Use freqency feature set to classify epileptiform event
    maxFrequency = double(max(spikeFrequency(:,2)));    %Calculate Maximum frequency during the SLE
    thresholdTonicFrequency = maxFrequency/3;

    if thresholdTonicFrequency < 1
        thresholdTonicFrequency = 1;    %set floor frequency to >1 hz (justify with data analysis)
    end
            
    indexTonic = spikeFrequency(:,2) >= thresholdTonicFrequency; %Use Michael's threshold to seperate frequency feature set into two populations, high and low; with floor frequency at 1 Hz
    spikeFrequency(:,3) = indexTonic;    %store Boolean index 

    %locate start of Tonic phase | Contingous segments above threshold    
    for j = 2: numel (indexTonic) %slide along the SLE; ignore the first spike, which is the sentinel spike
        if j+1 > numel(indexTonic)  %If you scan through the entire SLE and don't find a tonic phase, classify event as a IIE            
            j = find(indexTonic(2:end),1,'first');  %Take the 1st sec where frequency is 'high' as the onset if back-to-back high frequency are not found
                if isempty(j)
                    j = 1;  %honory position, just to push the function through
                end                
            startTonic(1,1) = spikeFrequency(j+1);  %j+1 because you started the find function (line above) starting from position "2:end"
            startTonic(1,2) = j+1;    %store the index
            endTonic(1,1) = spikeFrequency(numel(indexTonic));    %if no tonic period is found; just state whole ictal event as a tonic phase
            endTonic(1,2) = numel(indexTonic);  %store the index
            classification = 0; % 0 = no tonic phase, 1 = tonic-clonic SLE, 2 - tonic-only        
        else                        
            if indexTonic(j) > 0 && indexTonic(j+1) > 0 %If you locate two continuous segments with high frequency, mark the first segment as start of tonic phase                        
                startTonic(1,1) = spikeFrequency(j);  %store the onset time
                startTonic(1,2) = j;    %store the index
                classification = 1;   %1 = tonic-clonic SLE;   2 = tonic-only
                while ~and(indexTonic(j) == 0, indexTonic(j+1) == 0) & spikeFrequency(j,2) ~= 0; %Locate the offset time as the first point as either two continueous segments with low frequency or zero frequency, which ever comes first
                    j = j+1;    %keep sliding along the SLE until the statement above is false.
                    if j+1 > numel(indexTonic)  %If you slide all the way to the end and still can't find tonic phase offset, 
                        j = numel(indexTonic)+1;  %take the last point as the offset of the tonic phase - this means there is no clonic phase; add 1 because note you will remove it in line 33
                        classification = 2;   %1 = SLE;   2 = Tonic-only     
                        break
                    end                                    
                end            
                endTonic(1,1) = spikeFrequency(j-1); %take the point before the frequency drops to be the end of the tonic-phase
                endTonic(1,2) = j-1;  %store the index   
                %This is to confirm the tonic phase is >3s
                if (endTonic(1)-startTonic(1)) < (3*frequency)  %There needs to be at least 3 seconds of tonic phase
                    continue
                else                                  
                    break
                end                
            end
        end        
    end          

    %Characterize the epileptiform event | Label all the features of the SLE
    startTime = spikeFrequency(1,1)/frequency;  %secs
    endTime = spikeFrequency(end,1)/frequency;  %secs
    startTonicTime = startTonic(1)/frequency;  %secs
    endTonicTime = endTonic(1)/frequency;  %secs

    preictalPhaseFrequency = mean(spikeFrequency (1:startTonic(1,2),2));
    tonicPhaseFrequency = mean(spikeFrequency (startTonic(1,2):endTonic(1,2),2));
    minTonicPhaseFrequency = min(spikeFrequency (startTonic(1,2):endTonic(1,2),2));
    clonicPhaseFrequency = mean(spikeFrequency (endTonic(1,2):end,2));
    
%     %Tonic Phase must be >3 s; the onset trigger alone cannot be the tonic phase
%     if (endTonicTime-startTonicTime) <=3
%         classification = 0;
%     end
        
%     preictalPhaseDuration = startTonicTime - startTime;
%     tonicPhaseDuration = endTonicTime  - startTonicTime;
%     clonicPhaseDuration = endTime - endTonicTime; 
else    %this event is a IIS mascurading as a IIE
    startTime = spikeFrequency(1,1)/frequency;  %secs
    endTime = spikeFrequency(end,1)/frequency;  %secs
    startTonicTime = startTime;  %secs; just to push it through the function
    endTonicTime = endTime;  %secs; just to push it through the function
    
    classification = -1;    %-1 = IIS (<3s duration)  
    %fprintf(2,'\nWarning: The event at %d-%d s is not a IIE (duration <3 s), reclassify as a IIS.\n', startTime, endTime)  %Warning presented when a reclassification is made
    
    preictalPhaseFrequency = 0;
    tonicPhaseFrequency = 0;
    minTonicPhaseFrequency = 0;
    clonicPhaseFrequency = 0;

%     preictalPhaseDuration = startTonicTime - startTime;
%     tonicPhaseDuration = endTonicTime  - startTonicTime;
%     clonicPhaseDuration = endTime - endTonicTime; 
end
    
%Store event's characteristics for output
eventPhases(1) = classification;
eventPhases(2) = preictalPhaseFrequency;    %avg
eventPhases(3) = tonicPhaseFrequency;   %avg
eventPhases(4) = clonicPhaseFrequency;  %avg
eventPhases(5) = minTonicPhaseFrequency;    %the min value
eventPhases(6) = startTonicTime;
eventPhases(7) = endTonicTime;



