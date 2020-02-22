%% CWT script
%Instruction: Run analysisAlgorithm.m, then pause the function
%"dominantFrequency.m" at line 241

i=7

    %Event Vector
    eventVector = epileptiformEvent{i, 1};
    
    %Time Vector
    timeVector = (0:(length(eventVector)- 1))/frequency;
    timeVector = timeVector';    
    
    
    %Deciminated Event 
    eventVector = decimate(eventVector, 100);
    timeVector = decimate(timeVector,100);
    frequency_decimated = frequency/100;
    
    %Time Vector (decimated)
    timeVector = (0:(length(eventVector)- 1))/frequency_decimated;
    timeVector = timeVector';    
    
    
    %Frequency content of epileptiform event 
%     [s,f,t,p] = spectrogram (eventVector,round(windowSize*frequency),round(windowOverlap*frequency), 2.^nextpow2(windowSize*frequency), frequency, 'yaxis');
    [wt, f] = cwt(eventVector, 'amor', frequency);
    p = abs(wt);
    contourf(timeVector, f,p)
    
    
    %calculate decimated event vector
    [wt, f] = cwt(eventVector, 'amor', frequency_decimated);
    p = abs(wt);
    contourf(timeVector, f,p)
    
    
    
    figure()
    surface(timeVector, f,p)
    shading flat
    
    
    %Dominant Frequency at each time point 
    [maxS, idx] = max(p);        
    maxFreq = f(idx);   %finding the frequency with the maximum PSD
    indexInvalid = find (maxFreq > fc);  %Ignore all frequency above fc (frequency cutoff)
    maxFreq(indexInvalid) = 0;
    
    %store the max frequency content of each event 
    epileptiformEvent{i, 2} = t;
    epileptiformEvent{i, 3} = maxFreq;    

    %decipher
    [label, classification] = decipher (events,i);
    
    %Plot Figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('%s Event #%d', label, i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));

    subplot (3,1,1)
    plot (timeVector(round(windowOverlap*frequency):round(t(end)*frequency)), eventVector(round(windowOverlap*frequency):round(t(end)*frequency)))  %Plot ictal event
    hold on
    plot (timeVector(round(windowSize*frequency)), eventVector(round(windowSize*frequency)), 'ro', 'color', 'black', 'MarkerFaceColor', 'green')    %SLE onset
    plot (timeVector(round(numel(eventVector)-(windowSize*frequency))), eventVector(round(numel(eventVector)-(windowSize*frequency))), 'ro', 'color', 'black', 'MarkerFaceColor', 'red')    %SLE offset
    title (sprintf('LFP Bandpass Filtered (%s), %s Event #%d   |   Treatment Group:%d', filter, label, i, events(i,4)))
    xlabel('Time (sec)')
    ylabel('Voltage (mV)')
    axis tight   
    
    subplot (3,1,2)
    imagesc(t,f,10*log10(p))
    c = colorbar;
    c.Label.String = 'Power Spectral Density (dB/Hz)';  %Originally I wrote 'Power (dB)', but I think I've been calculating the PSD
    ylim([0 100])
    title (sprintf('Frequency Content (PSD) of %s Event #%d. Michaels Algorithm detected: %s', label, i, classification))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')         
    
    subplot (3,1,3)
    plot(t,maxFreq) 
    title (sprintf('Dominant Frequency over duration of %s Event #%d', label, i))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    axis tight
    ylim  ([0 70])   
    