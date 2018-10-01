function [figHandle] = plotEvent(timeSeries, eventTime,frequency)
%plotEvent plots vectors of event vectors are detected by your algorithm
%   plotEvent function can also plot addition markers and features about
%   the event that were detected by the detection algorithm. There is also
%   an option to plot baseline recording preceding and following the event
%   to provide 'context' of the baseline around the event. timeSeries is
%   the LFP signal to be plotted. The background bvector is the event
%   vector with 'context'.

%% Set variables to default values, if not specified
if nargin < 2
    frequency = 10000;    
end


  %% Plotting out detected Events with context     
    %make Event Vector
    onsetTime = (round(eventTime(1,1)*frequency));
    offsetTime = (round(eventTime(1,2)*frequency));
    eventVector = (onsetTime:offsetTime);  %SLE Vector    

    %make Background Vector
    if (onsetTime >= (context*frequency)+1 && (offsetTime+(context*frequency))<numel(timeSeries))
        backgroundVector = int64(onsetTime-50000:offsetTime+50000);   %Background Vector
    elseif (onsetTime < (context*frequency)+1)  %plus 1, because if the onsetTime happen to be exactly the same as, i.e. 5s, which is the preceding context, the starting index for the background vector will be 0 and cause an error
        backgroundVector = int64(1:offsetTime+(context*frequency));
    elseif ((offsetTime+(context*frequency))>numel(timeSeries))
        backgroundVector = int64(onsetTime-(context*frequency):numel(timeSeries));
    end

    %Plot figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('%s Event #%d', finalTitle, i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));  

    centerLFP = (timeSeries(backgroundVector(1)));  %center the LFP 
    plot (t(backgroundVector),timeSeries(backgroundVector)-centerLFP ) %background
    hold on    
    plot (t(eventVector),timeSeries(eventVector)-centerLFP)     %SLE
    plot (t(onsetTime), timeSeries(onsetTime)-centerLFP , 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
    plot (t(offsetTime), timeSeries(offsetTime)-centerLFP , 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
end


    

