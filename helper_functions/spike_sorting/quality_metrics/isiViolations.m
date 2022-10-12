function [fpRate, num_violations] = isiViolations(spike_train, min_time, max_time, isi_threshold, min_isi)
%Calculate ISI violations for a spike train.
%Based on metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
%modified by Dan Denman from cortex-lab/sortingQuality GitHub by Nick Steinmetz
%Inputs:
%-------
%spike_train : array of spike times
%min_time : minimum time for potential spikes
%max_time : maximum time for potential spikes
%isi_threshold : threshold for isi violation
%min_isi : threshold for duplicate spikes
%Outputs:
%--------
%fpRate : rate of contaminating spikes as a fraction of overall rate
%A perfect unit has a fpRate = 0
%A unit with some contamination has a fpRate < 0.5
%A unit with lots of contamination has a fpRate > 1.0
% num_violations : total number of violations

if nargin < 4
    isi_threshold = 0.0015;
end
if nargin < 5
    min_isi = 0;
end

duplicate_spikes = [0; diff(spike_train) <= min_isi];

spike_train = spike_train(~duplicate_spikes);
isis = diff(spike_train);

num_spikes = numel(spike_train);
num_violations = sum(isis < isi_threshold);
violation_time = 2*num_spikes*(isi_threshold - min_isi);
total_rate = firingRate(spike_train, min_time, max_time);
violation_rate = num_violations/violation_time;
fpRate = violation_rate/total_rate;

end