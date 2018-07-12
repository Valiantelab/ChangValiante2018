function [ L, LperSecond] = findLineLengths(samplingInterval, startTimes, endTimes, x, y)
    %findLineLengths: This function calculates the line length for a function
    %between specified time points
        %  x is the independent variable vector (i.e. time)
        %  y is the dependent variable vector you eant to find the line length of
        %  between specified start and end points
        %  Specify all the start times in the column vector start times, and all
        %  the end times in the column vector endTimes. Make sure these vectors
        %  have the same length for the function to work
        %  samplingInterval converts the time values from seconds to an index for
        %  the function to work
        %  The function returns a column vector L of line lengths of the function
        %  in between the specified time segments

    if length(startTimes) == length(endTimes)
        %Convert times from seconds to the index using the sampling interval. Note: seconds = (index - 1) / (sampling interval)
        startIndices = samplingInterval * startTimes + 1;
        endIndices = samplingInterval * endTimes + 1;
        
        %Initialize cells to store the dx, dy, dL values for each time segment. The number of cell columns corresponds to the number of start and end points provided by the user
        dx = cell(1,length(startIndices));
        dy = cell(1,length(startIndices));
        dL = cell(1,length(startIndices));
        
        L = NaN(length(startIndices), 1); %store L values for each time segment in a column vector
        duration = NaN(length(startIndices), 1); %store duration values for each time segment in a column vector
        
        %Calculate, dx, dy, dL values for each cell column. Store final L values in a column vector 
        for i = 1:length(startIndices)
            dx{:, i} = diff(x(startIndices(i) : endIndices(i)));
            dy{:, i} = diff(y(startIndices(i) : endIndices(i)));
            dL{:, i} = sqrt(dx{:, i}.^2 + dy{:, i}.^2);
            
            L(i) = sum(dL{:, i});
            duration(i) = endTimes(i) - startTimes(i);
        end
        
        LperSecond = L./duration;
    end
end