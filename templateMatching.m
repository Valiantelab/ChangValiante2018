% A template is given
% template = randn(100,1);
template = LFP(170900:189700);

figure;
plot(template)
axis tight
title ('Template')

% Create a matched filter based on the template
b = template;

% For testing the matched filter, create a random signal which
% contains a match for the template at some time index
% x = [randn(200,1); template(:); randn(300,1)];
x = LFP(1:300000);
n = 1:length(x);

figure;
plot(x)
axis tight
title ('Time Series Signal')

% Process the signal with the matched filter
y = filter(b,1,x);

% Set a detection threshold (exmaple used is 90% of template)
thresh = .9;

% Compute normalizing factor
u = template.'*template;

% Find matches
%matches = n(y>thresh*u);
matches = n(y>thresh*u);


% Plot the results
plot(n,y,'b', n(matches), y(matches), 'ro');

% Print the results to the console
% display(matches);