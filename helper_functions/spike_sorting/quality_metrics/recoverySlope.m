function recovery_slope = recoverySlope(waveform, timestamps, window)

%     """
%     Spike recovery slope (after repolarization)
%     Inputs:
%     ------
%     waveform : numpy.ndarray (N samples)
%     timestamps : numpy.ndarray (N samples)
%     window : int
%         Window (in samples) for linear regression
%     Outputs:
%     --------
%     recovery_slope : slope of recovery period (V / s)
%     """

if nargin < 3
    window = 20;
end

[~, max_point] = max(abs(waveform));

waveform = -1.*waveform * (sign(waveform(max_point)));

[~, peak_idx] = max(waveform(max_point:end));
peak_idx = peak_idx + max_point;

recovery_slope = timestamps(peak_idx:peak_idx+window) \ waveform(peak_idx:peak_idx+window);
recovery_slope = recovery_slope .* 1e-6;

end