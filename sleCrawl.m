%% Determine exact onset and offset times | Power Feature

%Need to write a function to determine when power increases, that's the
%point of seizure onset

% for i = 1:size(epileptiformLocation,1)
%     epileptiformLocation(i,1)


%% Scanning Low-Pass Filtered Power signal for more accurate onset/offset times

%testing for 25 hz low pass filter
for i = 1:size(SLE,1)
    
    %Rough SLE onset and offset times,  
    onsetSLE = int64((SLE(i,1)*10000));
    offsetSLE = int64((SLE(i,2))*10000);

    %SLE "context" (pre/post- baseline)
    onsetBaselineStart = (onsetSLE-10000);
    onsetBaselineEnd = (onsetSLE+5000);
    offsetBaselineStart = (offsetSLE-5000);
    offsetBaselineEnd = (offsetSLE+10000);

    %Range of LFP to search
    onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
    offsetContext = int64(offsetBaselineStart:offsetBaselineEnd);
    
    %Initial spike should be 1/3 the maximum prominance to be considered the onset 
    prominence = max(powerFeatureLowPassFiltered25(onsetContext))/3;
    
    %Locating potential onset points  
    [onset_pks, onset_locs, width_spike] = findpeaks(powerFeatureLowPassFiltered25(onsetContext), 'MinPeakProminence', prominence);
       
       
    %First spike is the onset
    SLEonset_final(i,1) = t(onsetContext(onset_locs(1)));
    
%     %Correction Factor
%     correctionFactor(i,1) = (width_spike(1)/2)/10000;   
%    
%     %Locating accurate onset time 
%     SLEonset_final(i,1) = SLEonset(i,1) - correctionFactor(i,1);
     
    
    %Test plot onset
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE onset #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    subplot (2,1,1)
    plot(t(onsetContext),LFP_normalized(onsetContext))
    hold on
    plot(t(onsetSLE), LFP_normalized(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
    plot(SLEonset(i,1), LFP_normalized(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14) 
    plot(t(onsetContext(onset_locs)), LFP_normalized(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('LFP normalized');
    ylabel ('mV');
    xlabel ('Time (sec)');
    
    subplot (2,1,2)
    plot(t(onsetContext), powerFeatureLowPassFiltered25(onsetContext))
    hold on
    plot(t(onsetSLE), powerFeatureLowPassFiltered25(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
    plot(SLEonset(i,1), powerFeatureLowPassFiltered25(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)
    plot(t(onsetContext(onset_locs)), powerFeatureLowPassFiltered25(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('Power, Low Pass Filtered (2 Hz)');
    ylabel ('mV');
    xlabel ('Time (sec)');
    
    
    %Test plot offset
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE onset #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    subplot (2,1,1)
    plot(t(onsetContext),LFP_normalized(onsetContext))
    hold on
    plot(t(onsetSLE), LFP_normalized(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
    plot(SLEonset(i,1), LFP_normalized(onsetContext(indexOnset)), 'o', 'color', 'black', 'MarkerSize', 14) 
    plot(t(onsetContext(onset_locs)), LFP_normalized(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('LFP normalized');
    ylabel ('mV');
    xlabel ('Time (sec)');
    
    subplot (2,1,2)
    plot(t(onsetContext), powerFeatureLowPassFiltered(onsetContext))
    hold on
    plot(t(onsetSLE), powerFeatureLowPassFiltered(onsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
    plot(SLEonset(i,1), powerFeatureLowPassFiltered(onsetContext(indexOnset)), 'o', 'color', 'black', 'MarkerSize', 14)
    plot(t(onsetContext(onset_locs)), powerFeatureLowPassFiltered(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
    %Labels
    title ('Power, Low Pass Filtered (2 Hz)');
    ylabel ('mV');
    xlabel ('Time (sec)');
end



    %Locating accurate offset time | peak finder | you also go to make sure it's not light-triggered  

for i = 1:size(SLE,1);
    
    
    %Rough SLE onset and offset times,  
    onsetSLE = int64((SLE(i,1)*10000));
    offsetSLE = int64((SLE(i,2))*10000);

    %SLE "context" (pre/post- baseline)
    onsetBaselineStart = (onsetSLE-10000);
    onsetBaselineEnd = (onsetSLE+5000);
    offsetBaselineStart = (offsetSLE-5000);
    offsetBaselineEnd = (offsetSLE+10000);

    %Range of LFP to search
    onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
    offsetContext = int64(offsetBaselineStart:offsetBaselineEnd);
    
    
    %make sure it's not light triggered
    %[Place Holder]
    
    %find the last spike    
    meanOffsetBaseline = mean(powerFeatureLowPassFiltered(offsetContext));
    OffsetLocation = powerFeatureLowPassFiltered(offsetContext) > meanOffsetBaseline/2;
    
    %Locating accurate offset time 
    indexOffset = find(OffsetLocation, 1, 'last');       
    SLEoffset(i,1) = t(offsetContext(indexOffset));

    
    %test plot offset
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE offset #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    subplot (2,1,1)
    plot(t(offsetContext),AbsLFP_normalizedFiltered(offsetContext))
    hold on
    plot(t(offsetSLE), AbsLFP_normalizedFiltered(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)
    plot(SLEoffset(i,1), AbsLFP_normalizedFiltered(offsetContext(indexOffset)), 'o', 'color', 'red', 'MarkerSize', 14)
    subplot (2,1,2)
    plot(t(offsetContext), powerFeatureLowPassFiltered(offsetContext))
    hold on
    plot(t(offsetSLE), powerFeatureLowPassFiltered(offsetSLE), 'x', 'color', 'red', 'MarkerSize', 12)
    plot(SLEoffset(i,1), powerFeatureLowPassFiltered(offsetContext(indexOffset)), 'o', 'color', 'red', 'MarkerSize', 14)

end


SLE_event = [SLEonset, SLEoffset];  %the detected events


[pks, locs, w] = findpeaks(data2);
max(pks)
numel(pks)

figure
plot(data2/10000)
hold on
plot (data2(locs), 'x', 'color', 'r')

[pks, locs] = findpeaks (data2);


%% Plot detected SLEs 

data1 = LFP_normalized;
data1 = diffLFP_normalized;
data2 = AbsLFP_normalizedFiltered;
 
data1 = DiffLFP_normalizedFiltered;
data2 = LFP_normalizedLowPassFiltered;

data2 = (DiffLFP_normalizedFiltered).^2;


for i = 1:size(SLE,1)
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('SLE #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
   
    time1 = single(SLE(i,1)*10000);
    time2= single(SLE(i,2)*10000);
    rangeSLE = (time1:time2);
    rangeOverview = (time1-50000:time2+50000);
    
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

