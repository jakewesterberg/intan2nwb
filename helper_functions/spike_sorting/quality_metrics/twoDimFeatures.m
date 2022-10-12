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
sites_to_sample = sites_to_sample((sites_to_sample>0) .* (sites_to_sample<waveform.shape[0])]

    wv = waveform[sites_to_sample, :]

    trough_idx = np.argmin(wv, 1)
    trough_amplitude = np.min(wv, 1)

    peak_idx = np.argmax(wv, 1)
    peak_amplitude = np.max(wv, 1)

    duration = np.abs(timestamps[peak_idx] - timestamps[trough_idx])

    overall_amplitude = peak_amplitude - trough_amplitude
    amplitude = np.max(overall_amplitude)
    max_chan = np.argmax(overall_amplitude)

    points_above_thresh = np.where(overall_amplitude > (amplitude * spread_threshold))[0]
    
    if len(points_above_thresh) > 1:
        points_above_thresh = points_above_thresh[isnot_outlier(points_above_thresh)]

    spread = len(points_above_thresh) * site_spacing * 1e6

    channels = sites_to_sample - peak_channel
    channels = channels[points_above_thresh]

    trough_times = timestamps[trough_idx] - timestamps[trough_idx[max_chan]]
    trough_times = trough_times[points_above_thresh]

    velocity_above, velocity_below = get_velocity(channels, trough_times, site_spacing)