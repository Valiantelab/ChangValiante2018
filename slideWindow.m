function [idx,result] = slideWindow(func, data, windowLength, overlap)

    if nargin == 3
        overlap = 0;
    end 
    
    windowEnd = windowLength;
    idx = [];
    result = [];
    while (windowEnd < length(data))
        idx = [idx  windowEnd];
        currentWindow = data(windowEnd-windowLength+1:windowEnd);
        
        result = [result func(currentWindow)];
        
        
        
        windowEnd = windowEnd + windowLength-overlap;
    end
    
    
    
end