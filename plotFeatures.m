figure

subplot (5,1,1)
plot(LFP_normalizedFiltered(onsetContext))
title('LFP filtered')

subplot (5,1,2)
plot(onsetContext, abs(LFP_normalizedFiltered(onsetContext)))
hold on
spikeIndex = find(onsetBaselineStart<=locs_spike & onsetBaselineEnd >=locs_spike);
for i = 1:numel(spikeIndex)
    plot(locs_spike(spikeIndex(i)), abs(LFP_normalizedFiltered((locs_spike(spikeIndex(i))))),'*r')
end
title('Absolute values of LFP filtered')

subplot (5,1,3)
plot (DiffLFP_normalizedFiltered(onsetContext))
title('Deriveative of LFP filtered')

subplot (5,1,4)
plot (powerFeature(onsetContext))
title('Power of LFP derivative ')

subplot (5,1,5)
plot (powerFeatureLowPassFiltered25(onsetContext))
title('power of LFp derivative, low pass filtered 25 Hz')
