function [average, sigma, Q] = quartilesStat (data)
%% statistics


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

