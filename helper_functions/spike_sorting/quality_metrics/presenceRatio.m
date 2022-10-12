function ratios = presenceRatio(spike_train)

h = histcounts(spike_train, linspace(min(spik_train), max(spike_train), 100));
ratios = sum(h>0) / 100;

end
