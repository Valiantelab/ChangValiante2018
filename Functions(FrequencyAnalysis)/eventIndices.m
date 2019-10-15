function [IndicesEvent,indicesBackground] = eventIndices(timeSeriesLFP, eventTime, contextDuration, frequency)
%indices.m function helps to create the indices of event vectors and
%background vectors for analysis.
%   This function is helpful because it helps to create the background
%   vectors of the first and last event, in the rare case that there is not
%   enough background "contextDuration" baseline preceding the first event or
%   subsequent to the last event.

if nargin <3
    contextDuration = 5;  %secs    
    frequency = 10000;  %Hz
end

%%Locating the indices of events to analyze
%Event Indices
onsetTime = (round(eventTime(1,1)*frequency));
offsetTime = (round(eventTime(1,2)*frequency));
IndicesEvent = (onsetTime:offsetTime);  %SLE Vector    

%Background Indices
if (onsetTime >= (contextDuration*frequency)+1 && (offsetTime+(contextDuration*frequency))<numel(timeSeriesLFP))
    indicesBackground = int64(onsetTime-(contextDuration*frequency):offsetTime+(contextDuration*frequency));   %Background Vector
elseif (onsetTime < (contextDuration*frequency)+1)  %plus 1, because if the onsetTime happen to be exactly the same as, i.e. 5s, which is the preceding contextDuration, the starting index for the background vector will be 0 and cause an error
    indicesBackground = int64(1:offsetTime+(contextDuration*frequency));
elseif ((offsetTime+(contextDuration*frequency))>numel(timeSeriesLFP))
    indicesBackground = int64(onsetTime-(contextDuration*frequency):numel(timeSeriesLFP));
end

end

