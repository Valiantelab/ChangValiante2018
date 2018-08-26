%% Plotting out the detected SLEs | To figure out how off you are
%define variables
data1 = LFP_normalized;
data2 = powerFeature;

for i = 1:size(SLE_final,1)
%     figure;
%     set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('SLE #%d', i)); %select the name you want
%     set(gcf, 'Position', get(0, 'Screensize'));   
%    
    time1 = single(SLE_final(i,1)*10000);
    time2= single(SLE_final(i,2)*10000);
    rangeSLE = (time1:time2);
    SLE_Vector{i,1} = rangeSLE;   %store the datarange for SLEs

    if time1>50001
        rangeOverview = (time1-50000:time2+50000);
    else
        rangeOverview = (1:time2+50000);
    end
    
    SLE_Vector{i,2} = rangeOverview;   %store the SLE vectors
 
end

    
    subplot (2,1,1)
    %overview
    plot (t(rangeOverview),data1(rangeOverview))
    hold on
    %SLE
    plot (t(rangeSLE),data1(rangeSLE))
    %onset/offset markers
    plot (t(time1), data1(time1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset
    plot (t(time2), data1(time2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset
    %Labels
    title ('Derivative of Filtered LFP');
    ylabel ('mV');
    xlabel ('Time (sec)');


    subplot (2,1,2)
    %overview
    plot (t(rangeOverview),data2(rangeOverview))
    hold on
    %SLE
    plot (t(rangeSLE),data2(rangeSLE))
    %onset/offset markers
    plot (t(time1), data2(time1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset
    plot (t(time2), data2(time2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset
    %Labels
    title ('Absolute Derivative of Filtered LFP');
    ylabel ('mV');
    xlabel ('Time (sec)');
end