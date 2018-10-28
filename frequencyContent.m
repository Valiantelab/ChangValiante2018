vector = LFP_filtered((events(end,1)*frequency):(events(end,2)*frequency));
timeVector = (0:(length(vector)-1))/frequency;
timeVector = timeVector';
% figure
% plot (vector)

% vector_d = decimate (vector, 100);
% figure
% [s,f,t] = spectrogram (vector_d, 100, 80, [], 100, 'yaxis');

[s,f,t] = spectrogram (vector, 100000, 90000, [], 10000, 'yaxis');

% figure;
% sdB = 10*log10(abs(s).^2);
% contourf(t,f,sdB)
% ylim ([0 30])
% mesh(t,f,sdB)
% c = colorbar;

figure;
subplot (3,1,1)
plot (timeVector, vector)
title (sprintf('LFP Bandpass Filtered (0-100 Hz), SLE #%d', i))
xlabel('Time (sec)')
ylabel('Voltage (mV)')
axis tight

subplot (3,1,2)
contour(t,f,abs(s).^2)
c = colorbar;
c.Label.String = 'Power (Watt)'
ylim([0 20])
i=1
title (sprintf('Frequency Content of SLE #%d', i))
ylabel('Frequency (Hz)')
xlabel('Time (sec)')

%%max frequency
subplot (3,1,3)
[maxS, idx] = max(abs(s));
maxFreq = f(idx);
plot(t,maxFreq)
title (sprintf('Dominant Frequency over duration of SLE #%d', i))
ylabel('Frequency (Hz)')
xlabel('Time (sec)')
axis tight