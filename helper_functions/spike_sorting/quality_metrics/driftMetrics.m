function driftMetrics(spike_times, spike_clusters, spike_templates, total_units, ...
    pc_features,pc_feature_ind,interval_length,min_spikes_per_interval)

%         """
%         Helper to calculate drift for one cluster
%         Args:
%             cluster_id:
%         Returns:
%             max_drift, cumulative_drift
%         """

        in_cluster = spike_clusters == cluster_id
        times_for_cluster = spike_times[in_cluster]
        depths_for_cluster = depths[in_cluster]

        median_depths = []

        for t1, t2 in zip(interval_starts, interval_ends):

            in_range = (times_for_cluster > t1) * (times_for_cluster < t2)

            if np.sum(in_range) >= min_spikes_per_interval:
                median_depths.append(np.median(depths_for_cluster[in_range]))
            else:
                median_depths.append(np.nan)

        median_depths = np.array(median_depths)
        max_drift = np.around(np.nanmax(median_depths) - np.nanmin(median_depths), 2)
        cumulative_drift = np.around(np.nansum(np.abs(np.diff(median_depths))), 2)
        return max_drift, cumulative_drift


    max_drifts = []
    cumulative_drifts = []

    depths = get_spike_depths(spike_templates, pc_features, pc_feature_ind)

    interval_starts = np.arange(np.min(spike_times), np.max(spike_times), interval_length)
    interval_ends = interval_starts + interval_length

    cluster_ids = np.unique(spike_clusters)

    if do_parallel:
        from joblib import Parallel, delayed
        meas = Parallel(n_jobs=-1, verbose=2)(delayed(calc_one_cluster)(cluster_id)
                                              for cluster_id in cluster_ids)
    else:
        meas = [calc_one_cluster(cluster_id) for cluster_id in cluster_ids]

    for max_drift, cumulative_drift in meas:
        max_drifts.append(max_drift)
        cumulative_drifts.append(cumulative_drift)
    return np.array(max_drifts), np.array(cumulative_drifts)