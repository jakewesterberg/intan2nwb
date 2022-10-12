function wh = waveHalfwidth(waveform, timestamps)
% 
% """
% Spike width (in seconds) at half max amplitude
% Inputs:
% ------
% waveform : numpy.ndarray (N samples)
% timestamps : numpy.ndarray (N samples)
% Outputs:
% --------
% halfwidth : waveform halfwidth in milliseconds
% """

[~,trough_idx] = min(waveform);
[~,peak_idx] = max(waveform);

try
    if waveform(peak_idx) > abs(waveform(trough_idx))
        threshold = waveform(peak_idx) * 0.5;
        thresh_crossing_1 = find(waveform(1:peak_idx) > threshold,1);
        thresh_crossing_2 = find(waveform(peak_idx:end) < threshold,1) + peak_idx;
    else
        threshold = waveform(trough_idx) * 0.5;
        thresh_crossing_1 = find(waveform(1:trough_idx) < threshold,1);
        thresh_crossing_2 = find(waveform(trough_idx:end) > threshold,1) + trough_idx;
    end
    wh = timestamps(thresh_crossing_2) - timestamps(thresh_crossing_1);
catch
    wh = NaN;
end

wh = wh * 1e3;
end
