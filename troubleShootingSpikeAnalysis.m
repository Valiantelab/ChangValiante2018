%% Spiking Frequency Classifier
data1 = LFP_normalized; %Time series to be plotted 
lightpulse = LED > 1;

%for i = 1:size(putativeSLE,1)   
for i = 3
    %make SLE vector
    onsetTime = single(putativeSLE(i,1)*10000);
    offsetTime = single(putativeSLE(i,2)*10000);
    sleVector = (onsetTime:offsetTime);  %SLE Vector  
    SLE_vector{i} = sleVector;  %store SLE vector
    
    %Calculate the spiking rate for SLE
    windowSize = 1;  %seconds  
    
    sleDuration = round(numel(sleVector)/frequency);    %rounded to whole number
    clear spikeRateMinute
    for j = 1:sleDuration
        startWindow = onsetTime+((windowSize*frequency)*(j-1));
        EndWindow = onsetTime+((windowSize*frequency)*j);
        spikeRate = and(startWindow<=locs_spike_2nd, EndWindow >=locs_spike_2nd);
        spikeRateMinute(j,1) = startWindow;
        spikeRateMinute(j,2) = sum(spikeRate(:));   %number of spikes in the window
    end
    
    %average spiking rate of SLE
    putativeSLE (i,4) = mean(spikeRateMinute(:,2));
    
    %average spiking rate in 1st half, SLE
    putativeSLE (i,5) = mean(spikeRateMinute(int64(1:(sleDuration/2)),2));
    
    %average spiking rate in 2nd half, SLE
    putativeSLE (i,6) = mean(spikeRateMinute(int64((sleDuration/2):sleDuration),2));    
    
    %average power
    totalPower = sum(powerFeature(sleVector));
    putativeSLE (i,7) = totalPower /sleDuration;  
    
%     %Store SLE
%     SLE = putativeSLE(putativeSLE(:,4)>0, 1);
        
%     %calculate the average spiking rate at defined window
%     averageWindowSize = 5;  %seconds
%     x  = spikeRateMinute(:,2);
%     S  = numel(x);
%     xx = reshape(x(1:averageWindowSize - mod(S, averageWindowSize)), averageWindowSize, []);
%     y  = sum(xx, 1).' / averageWindowSize;
%         
    %make background vector
    if onsetTime >= 50001
        backgroundVector = (onsetTime-50000:offsetTime+50000);   %Background Vector
    else
        backgroundVector = (1:offsetTime+50000);
    end
    background_vector{i} = backgroundVector;  %store background vector

    %% plot vectors
    %if threshold_multiple(3) == 1   
    
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('Putative SLE #%d', i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));   
    
    plot (t(backgroundVector),data1(backgroundVector))
    hold on
    plot (t(sleVector),data1(sleVector))     %SLE
    plot (t(onsetTime), data1(onsetTime), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %onset marker
    plot (t(offsetTime), data1(offsetTime), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'red') %offset marker
    indexSpikes = and(onsetTime<locs_spike_2nd, offsetTime>locs_spike_2nd); %Locate spikes between the onset and offset  
    plot (t(locs_spike_2nd(indexSpikes)), (data1(locs_spike_2nd(indexSpikes))), 'x') %plot spikes (artifact removed)
    plot (t(backgroundVector),(lightpulse(backgroundVector)-1)/20, 'b') %plot LED   
    title (sprintf('Absolute values of Filtered LFP Recording, SLE #%d', i));
    ylabel ('mV');
    xlabel ('Time (sec)');   
    
    yyaxis right
    
    plot (spikeRateMinute(:,1)/frequency, spikeRateMinute(:,2), 'o', 'color', 'b')
    ylabel ('spike rate/second (Hz)');
    set(gca,'fontsize',16)
    %end
    
end