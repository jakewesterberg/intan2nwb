function fr = firingRate(spike_train, min_time, max_time)
% Calculate firing rate for a spike train.
% If no temporal bounds are specified, the first and last spike time are used.
% Inputs:
% -------
% spike_train : numpy.ndarray
% Array of spike times in seconds
% min_time : float
% Time of first possible spike (optional)
% max_time : float
% Time of last possible spike (optional)
% Outputs:
% --------
% fr : float
% Firing rate in Hz

if exist('min_time', 'var') & exist('max_time', 'var')
    duration = max_time - min_time;
else
    duration = max(spike_train) - min(spike_train);
end

fr = numel(spike_train) / duration;

end