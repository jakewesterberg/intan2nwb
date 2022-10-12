function [mean_waveforms, spike_count] = ...
    waveExtraction(raw_data, spike_times, spike_clusters, templates, ...
    channel_map, bit_volts, sample_rate, site_spacing, params)

%     """
%     Calculate mean waveforms for sorted units.
%     Inputs:
%     -------
%     raw_data : continuous data as numpy array (samples x channels)
%     spike_times : spike times (in samples)
%     spike_clusters : cluster IDs for each spike time []
%     clusterIDs : all unique cluster ids
%     cluster_quality : 'noise' or 'good'
%     sample_rate : Hz
%     site_spacing : m
%     Outputs:
%     -------
%     mean_waveforms : numpy array with dims :
%      - 1 : clusterID
%      - 2 : epochs
%      - 3 : mean (0) or std (1)
%      - 4 : channels
%      - 5 : samples
%     spike_count : numpy array with dims :
%      - 1 : clusterID
%      - 2 : epoch (last is entire dataset)
%     dimCoords : list of coordinates for each dimension
%     dimLabels : list of labels for each dimension
%     metrics : DataFrame with waveform metrics
%     Parameters:
%     ----------
%     samples_per_spike : number of samples in extracted spikes
%     pre_samples : number of samples prior to peak
%     num_epochs : number of epochs to calculate mean waveforms
%     spikes_per_epoch : max number of spikes to generate average for epoch
%     """
%
%     # #############################################


samples_per_spike       = 82; %'Number of samples to extract for each spike')
pre_samples             = 20; %'Number of samples between start of spike and the peak')
spikes_per_epoch        = 100; %'Max number of spikes per epoch')
upsampling_factor       = 200/82; %'Upsampling factor for calculating waveform metrics')
spread_threshold        = 0.12; %'Threshold for computing channel spread of 2D waveform')
site_range              = 16; %'Number of sites to use for 2D waveform metrics')

total_units             = numel(spike_clusters);
cluster_ids             = 1:total_units;

mean_waveforms          = zeros(total_units, 2, size(raw_data,1), samples_per_spike);
spike_count             = zeros(total_units, 1);

peak_channels           = squeeze(channel_map(np.argmax(max(templates) - min(templates),1)));

for cluster_idx = cluster_ids

    waveforms = nan(spikes_per_epoch, size(raw_data,1), samples_per_spike);

    np.random.shuffle(times_for_cluster)

    total_waveforms = min([times_for_cluster.size, spikes_per_epoch]);

    for wv_idx, peak_time in enumerate(times_for_cluster[:total_waveforms]):
        start = int(peak_time-pre_samples);
        end1 = start + samples_per_spike;
        rawWaveform = raw_data[start:end1, :].T
    end

    if rawWaveform.shape[1] == samples_per_spike
        waveforms[wv_idx, :, :] = rawWaveform * bit_volts;
    end

    with warnings.catch_warnings():

    warnings.simplefilter("ignore", category=RuntimeWarning)
    mean_waveforms[cluster_idx, epoch_idx,
        0, :, :] = np.nanmean(waveforms, 0)
    mean_waveforms[cluster_idx, epoch_idx,
        1, :, :] = np.nanstd(waveforms, 0)

    for channel in range(0, mean_waveforms.shape[3]):
        mean_waveforms[cluster_idx, 0, channel, :] = ...
            mean_waveforms[cluster_idx, 0, channel, :] - ...
            mean_waveforms[cluster_idx, 0, channel, 0];
    end

    spike_count[cluster_idx, epoch_idx] = total_waveforms;

end
end
