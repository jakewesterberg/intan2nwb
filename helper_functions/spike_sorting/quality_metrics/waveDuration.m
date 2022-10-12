function wd = waveDuration(waveform, timestamps)

%     """
%     Duration (in seconds) between peak and trough
%     Inputs:
%     ------
%     waveform : numpy.ndarray (N samples)
%     timestamps : numpy.ndarray (N samples)
%     Outputs:
%     --------
%     duration : waveform duration in milliseconds
%     """

[~,trough_idx] = min(waveform);
[~,peak_idx] = max(waveform);

if waveform(peak_idx) > abs(waveform(trough_idx))
    wd = timestamps(peak_idx:end);
    wd = wd(waveform(peak_idx:end) == min(waveform(peak_idx:end))) - timestamps(peak_idx);
else
    wd = timestamps(trough_idx:end);
    wd = wd(waveform(trough_idx:end) == max(waveform(trough_idx:end))) - timestamps(trough_idx);
end

wd = wd * 1e3;
end