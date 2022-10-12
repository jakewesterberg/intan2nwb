function twoDimFeatures(waveform, timestamps, peak_channel, spread_threshold = 0.12, site_range=16, site_spacing=10e-6):
    
%     """ 
%     Compute features of 2D waveform (channels x samples)
%     Inputs:
%     ------
%     waveform : numpy.ndarray (N channels x M samples)
%     timestamps : numpy.ndarray (M samples)
%     peak_channel : int
%     spread_threshold : float
%     site_range: int
%     site_spacing : float
%     Outputs:
%     --------
%     amplitude : uV
%     spread : um
%     velocity_above : s / m
%     velocity_below : s / m
%     """

    assert site_range % 2 == 0 # must be even

    sites_to_sample = np.arange(-site_range, site_range+1, 2) + peak_channel

    sites_to_sample = sites_to_sample[(sites_to_sample > 0) * (sites_to_sample < waveform.shape[0])]

    wv = waveform[sites_to_sample, :]

    #smoothed_waveform = np.zeros((wv.shape[0]-1,wv.shape[1]))
    #for i in range(wv.shape[0]-1):
    #    smoothed_waveform[i,:] = np.mean(wv[i:i+2,:],0)

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
 
    return amplitude, spread, velocity_above, velocity_below
