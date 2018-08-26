%% IIS: Determine exact onset and offset times | Power Feature
% Scanning Low-Pass Filtered Power signal for more accurate onset/offset times
for i = 1:size(IIS,1)
    
    %Rough IIS onset and offset times,  
    onsetIIS = int64((IIS(i,1)*10000));
    offsetIIS = int64((IIS(i,2))*10000);

    %IIS "context" (pre/post- baseline)
    onsetBaselineStart = (onsetIIS-10000);     
    onsetBaselineEnd = (onsetIIS+5000);
    offsetBaselineStart = (offsetIIS-5000);
    offsetBaselineEnd = (offsetIIS+50000);

    %Range of LFP to search
    onsetContext = int64(onsetBaselineStart:onsetBaselineEnd);
    offsetContext = int64(offsetBaselineStart:offsetBaselineEnd); 

    %Locating the onset time
    threshold = max(powerFeatureLowPassFiltered25(onsetContext))/4; %IIS onset where spike prominience > 1/4 the maximum amplitude
    [onset_pks, onset_locs, width_pks] = findpeaks(powerFeatureLowPassFiltered25(onsetContext), 'MinPeakHeight', threshold);     
    IISonset(i,1) = t(onsetContext(onset_locs(1))); %First spike is the onset   
   
    %Correction Factor
    correctionFactor = (width_pks(1)/2)/frequency;
    
    %Locating the corrected onset time
    IISonset_final(i,1) = IISonset(i,1)-correctionFactor;
    
    %Locating the offset time
    meanOffsetBaseline = mean(LFP_normalizedLowPassFiltered(offsetContext)); %IIS ends when signal returned to half the mean power of signal
    OffsetLocation = LFP_normalizedLowPassFiltered(offsetContext) > meanOffsetBaseline/1; 
    indexOffset = find(OffsetLocation, 1, 'last'); %Last point is the offset     
    IISoffset_final(i,1) = t(offsetContext(indexOffset)); 
    
%     %locating the offset time (test)
%     offsetThreshold = (LFP_normalizedLowPassFiltered(onsetContext(1))); %IIS ends when signal returned to half the mean power of signal
%     OffsetLocation = LFP_normalizedLowPassFiltered(offsetContext) > offsetThreshold; 
%     indexOffset_test = find(OffsetLocation, 1, 'last'); %Last point is the offset     
%     IISoffset_test(i,1) = t(offsetContext(indexOffset_test)); 
%     
    
    %plot onset
%     figure;
%     set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('IIS onset #%d', i)); %select the name you want
%     set(gcf, 'Position', get(0, 'Screensize'));   
%     subplot (2,1,1)
%     plot(t(onsetContext),LFP_normalized(onsetContext))
%     hold on
%     plot(t(onsetIIS), LFP_normalized(onsetIIS), 'x', 'color', 'red', 'MarkerSize', 12)  %initial (rough) detection
%     plot(IISonset_final(i,1), LFP_normalized(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)   %Detected onset point 
%     plot(t(onsetContext(onset_locs)), LFP_normalized(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     %Labels
%     title ('LFP normalized');
%     ylabel ('mV');
%     xlabel ('Time (sec)');
%     
%     subplot (2,1,2)
%     plot(t(onsetContext), LFP_normalizedLowPassFiltered(onsetContext))
%     hold on
%     plot(t(onsetIIS), LFP_normalizedLowPassFiltered(onsetIIS), 'x', 'color', 'red', 'MarkerSize', 12)     %initial (rough) detection
%     plot(IISonset_final(i,1), LFP_normalizedLowPassFiltered(onsetContext(onset_locs(1))), 'o', 'color', 'black', 'MarkerSize', 14)    %Detected onset point 
%     plot(t(onsetContext(onset_locs)), LFP_normalizedLowPassFiltered(onsetContext(onset_locs)), '*', 'color', 'green', 'MarkerSize', 14)    %Final (Refined) detection
%     %Labels
%     title ('Power, Low Pass Filtered (2 Hz)');
%     ylabel ('mV');
%     xlabel ('Time (sec)');
%     
%     
%     %plot offset
%     figure;
%     set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
%     set(gcf,'Name', sprintf ('IIS offset #%d', i)); %select the name you want
%     set(gcf, 'Position', get(0, 'Screensize'));   
%     subplot (2,1,1) %LFP overview
%     plot(t(offsetContext),LFP_normalized(offsetContext))
%     hold on
%     plot(t(offsetIIS), LFP_normalized(offsetIIS), 'x', 'color', 'red', 'MarkerSize', 12)     %initial detection as offset
%     plot(IISoffset_final(i,1), LFP_normalized(offsetContext(indexOffset)), 'o', 'color', 'green', 'MarkerSize', 14)
%     subplot (2,1,2) %feature extraction
%     plot(t(offsetContext), LFP_normalizedLowPassFiltered(offsetContext))
%     hold on
%     plot(t(offsetIIS), LFP_normalizedLowPassFiltered(offsetIIS), 'x', 'color', 'red', 'MarkerSize', 12)   %initial detection as offset
%     plot(IISoffset_final(i,1), LFP_normalizedLowPassFiltered(offsetContext(indexOffset)), 'o', 'color', 'green', 'MarkerSize', 14)
%     
end

IIS_final = [IISonset_final, IISoffset_final];  %final list of IISs, need to filter out artifacts

%% Plotting out the detected IISs | To figure out how off you are
%define variables
data1 = LFP_normalized;
data2 = LFP_normalizedLowPassFiltered;

for i = 1:size(IIS_final,1)
    figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('IIS #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
   
    time1 = int64(IIS_final(i,1)*10000);
    time2= int64(IIS_final(i,2)*10000);
    rangeIIS = int64((time1:time2));
    if time1>50001
        rangeOverview = int64((time1-50000:time2+50000));
    else
        rangeOverview = int64((1:time2+50000));
    end
    

    
    subplot (2,1,1)
    %overview
    plot (t(rangeOverview),data1(rangeOverview))
    hold on
    %IIS
    plot (t(rangeIIS),data1(rangeIIS))
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
    %IIS
    plot (t(rangeIIS),data2(rangeIIS))
    %onset/offset markers
    plot (t(time1), data2(time1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset
    plot (t(time2), data2(time2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset
    %Labels
    title ('Absolute Derivative of Filtered LFP');
    ylabel ('mV');
    xlabel ('Time (sec)');
end
