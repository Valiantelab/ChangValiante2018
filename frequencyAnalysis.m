function [figHandle] = frequencyAnalysis(timeSeries, events, time, i, frequency, troubleshooting)
%frequencyAnalysis produces thefrequency content of epileptiform events.
%   This function uses spectrogram to analyze the frequency content of the
%   epileptiform forms that are entered into the function. The window size
%   being analyzed is 10s, so the smallest event can only be 10 s (a SLE,
%   by definition). This allows the Rayleigh frequency to be 0.1 Hz.
%   Authors: Michael Chang, Liam Long, and Kramay Patel.

%Set default values, if not specified
if nargin <3
    i = 00;
    time = []
    frequency = 10000;  %Hz\
    troubleshooting = [];    
end

if troubleshooting 
  
    %Set variables    
    startLocation = int64(events(1)*frequency);
    endLocation = int64(events(2)*frequency);

    %Epileptiform events to analyze
    eventVector = timeSeries(startLocation:endLocation);  %event vector
    
    %Time vector
    if time 
        timeVector = time(startLocation:endLocation);   %make time Vector
    else
        timeVector = (0:(length(eventVector)-1))/frequency;  %make time vector
        timeVector = timeVector';
    end
    
    %Energy content of epileptiform event
    [s,f,t] = spectrogram (eventVector, 5*frequency, 4*frequency, [], frequency, 'yaxis');

    %Dominant Frequency at each time point
    [maxS, idx] = max(abs(s));
    maxFreq = f(idx);

    %decipher
    [label,classification] = decipher (events);
    
    %Plot Figures
    figHandle = figure;
    set(gcf,'NumberTitle','off', 'color', 'w'); %don't show the figure number
    set(gcf,'Name', sprintf ('%s Event #%d', label, i)); %select the name you want
    set(gcf, 'Position', get(0, 'Screensize'));
    
    subplot (3,1,1)
    plot (timeVector, eventVector)
    title (sprintf('LFP Bandpass Filtered (0-100 Hz), %s Event #%d', label, i))
    xlabel('Time (sec)')
    ylabel('Voltage (mV)')
    axis tight

    subplot (3,1,2)
    contour(t,f,abs(s).^2)
    c = colorbar;
    c.Label.String = 'Power (mV^2)';    %what is the unit really called? 
    ylim([0 100])
    set(gca, 'YScale', 'log')
    title (sprintf('Frequency Content of %s Event #%d. Michaels Algorithm detected: %s', label, i, classification))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')

    subplot (3,1,3)
    plot(t,maxFreq) 
    title (sprintf('Dominant Frequency over duration of %s Event #%d', label, i))
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    axis tight
    
end





