function [average, sigma, Q] = statistics (data)
% statistics function calculates the population statistics for the time
% series provided. Outputs the average value, standard deviation (sigma) of
% the values, and the quartile range for the time series.

% compute mean
average = mean(data);

% compute the standard deviation
sigma = std(data);

% compute the median
medianx = median(data);

% STEP 1 - rank the data
y = sort(data);

% compute 25th percentile (first quartile)
Q(1) = median(y(find(y<median(y))));

% compute 50th percentile (second quartile)
Q(2) = median(y);

% compute 75th percentile (third quartile)
Q(3) = median(y(find(y>median(y))));

end


