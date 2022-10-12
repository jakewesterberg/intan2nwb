function [amplitude, spread, velocity_above, velocity_below] = ...
    twoDimFeatures(waveform, timestamps, peak_channel, spread_threshold, site_range, site_spacing)

%     """
%     Compute features of 2D waveform (channels x samples)
%     Inputs:
%     ------
%     waveform : numpy.ndarray (N channels x M samples)
%     timestamps : numpy.ndarray (M samples)
%     peak_channel : int
%     spread_threshold : float
%     site_range: int and even
%     site_spacing : float
%     Outputs:
%     --------
%     amplitude : uV
%     spread : um
%     velocity_above : s / m
%     velocity_below : s / m
%     """

if nargin < 4
    spread_threshold = 0.12;
end
if nargin < 5
    site_range = 16;
end
if nargin < 6
    site_spacing = 10e-6;
end

sites_to_sample = (-1*site_range:2:site_range+1) + peak_channel;
sites_to_sample = sites_to_sample((sites_to_sample>0) .* (sites_to_sample<size(waveform,1)));

wv = waveform(sites_to_sample, :);

[trough_amplitude, trough_idx] = min(wv);
[peak_amplitude, peak_idx] = max(wv);

%wduration = abs(timestamps(peak_idx) - timestamps(trough_idx));

overall_amplitude = peak_amplitude - trough_amplitude;
[amplitude, max_chan] = max(overall_amplitude);

points_above_thresh = overall_amplitude > (amplitude * spread_threshold);

if length(points_above_thresh) > 1
    points_above_thresh = points_above_thresh(isNotOutlier(points_above_thresh));
end

spread = length(points_above_thresh) * site_spacing * 1e6;

channels = sites_to_sample - peak_channel;
channels = channels(points_above_thresh);

trough_times = timestamps(trough_idx) - timestamps(trough_idx(max_chan));
trough_times = trough_times(points_above_thresh);

[velocity_above, velocity_below] = getVelocity(channels, trough_times, site_spacing);

end