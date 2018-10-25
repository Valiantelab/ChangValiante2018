vector = LFP_filtered((events(end,1)*frequency):(events(end,2)*frequency));

figure
plot (vector_d)

vector_d = decimate (vector, 100);
figure
[s,f,t] = spectrogram (vector_d, 100, 80, [], 100, 'yaxis');

[s,f,t] = spectrogram (vector, 100000, 90000, [], 10000, 'yaxis');
figure;
sdB = 10*log10(abs(s).^2);
contourf(t,f,sdB)
ylim ([0 30])
mesh(t,f,sdB)
colorbar

contourf(t,f,abs(s))
lim([0 30])

%%max frequency
[maxS, idx] = max(abs(s));
maxFreq = f(idx);
plot(t,maxFreq)