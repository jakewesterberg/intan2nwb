function pcMetrics(spike_clusters,
                         spike_templates,
                         total_units,
                         pc_features,
                         pc_feature_ind,
                         num_channels_to_compare,
                         max_spikes_for_cluster,
                         max_spikes_for_nn,
                         n_neighbors,
                         do_parallel=True):
    """
    :param spike_clusters:
    :param total_units:
    :param pc_features:
    :param pc_feature_ind:
    :param num_channels_to_compare:
    :param max_spikes_for_cluster:
    :param max_spikes_for_nn:
    :param n_neighbors:
    :return:
    """

    assert (num_channels_to_compare % 2 == 1)
    half_spread = int((num_channels_to_compare - 1) / 2)

    cluster_ids = np.unique(spike_clusters)
    template_ids = np.unique(spike_templates)

    template_peak_channels = np.zeros((len(template_ids),), dtype='uint16')
    cluster_peak_channels = np.zeros((len(cluster_ids),), dtype='uint16')

    for idx, template_id in enumerate(template_ids):
        for_template = np.squeeze(spike_templates == template_id)
        pc_max = np.argmax(np.mean(pc_features[for_template, 0, :], 0))
        template_peak_channels[idx] = pc_feature_ind[template_id, pc_max]

    for idx, cluster_id in enumerate(cluster_ids):
        for_unit = np.squeeze(spike_clusters == cluster_id)
        templates_for_unit = np.unique(spike_templates[for_unit])
        template_positions = np.where(np.isin(template_ids, templates_for_unit))[0]
        cluster_peak_channels[idx] = np.median(template_peak_channels[template_positions])

    # Loop over clusters:
    if do_parallel:
        from joblib import Parallel, delayed
        meas = Parallel(n_jobs=-1, verbose=3)(  # -1 means use all cores
            delayed(calculate_pc_metrics_one_cluster)  # Function
            (cluster_peak_channels, idx, cluster_id, cluster_ids,
             half_spread, pc_features, pc_feature_ind,
             spike_clusters, spike_templates,
             max_spikes_for_cluster, max_spikes_for_nn, n_neighbors
             )
            for idx, cluster_id in enumerate(cluster_ids))  # Loop
    else:
        from tqdm import tqdm
        meas = []
        for idx, cluster_id in tqdm(enumerate(cluster_ids), total=cluster_ids.max(), desc='PC metrics'):  # Loop
            meas.append(calculate_pc_metrics_one_cluster(  # Function
                cluster_peak_channels, idx, cluster_id, cluster_ids,
                half_spread, pc_features, pc_feature_ind,
                spike_clusters, spike_templates,
                max_spikes_for_cluster, max_spikes_for_nn, n_neighbors))

    # Unpack:
    isolation_distances = []
    l_ratios = []
    d_primes = []
    nn_hit_rates = []
    nn_miss_rates = []
    for mea in meas:
        isolation_distance, d_prime, nn_miss_rate, nn_hit_rate, l_ratio = mea
        isolation_distances.append(isolation_distance)
        d_primes.append(d_prime)
        nn_miss_rates.append(nn_miss_rate)
        nn_hit_rates.append(nn_hit_rate)
        l_ratios.append(l_ratio)

    return (np.array(isolation_distances), np.array(l_ratios), np.array(d_primes),
            np.array(nn_hit_rates), np.array(nn_miss_rates))

def calculate_pc_metrics_one_cluster(cluster_peak_channels, idx, cluster_id,cluster_ids,
                                         half_spread, pc_features, pc_feature_ind,
                                         spike_clusters, spike_templates,
                                         max_spikes_for_cluster, max_spikes_for_nn, n_neighbors):

    peak_channel = cluster_peak_channels[idx]
    num_spikes_in_cluster = np.sum(spike_clusters == cluster_id)

    half_spread_down = peak_channel \
        if peak_channel < half_spread \
        else half_spread

    half_spread_up = np.max(pc_feature_ind) - peak_channel \
        if peak_channel + half_spread > np.max(pc_feature_ind) \
        else half_spread

    channels_to_use = np.arange(peak_channel - half_spread_down, peak_channel + half_spread_up + 1)
    units_in_range = cluster_ids[np.isin(cluster_peak_channels, channels_to_use)]

    spike_counts = np.zeros(units_in_range.shape)

    for idx2, cluster_id2 in enumerate(units_in_range):
        spike_counts[idx2] = np.sum(spike_clusters == cluster_id2)

    if num_spikes_in_cluster > max_spikes_for_cluster:
        relative_counts = spike_counts / num_spikes_in_cluster * max_spikes_for_cluster
    else:
        relative_counts = spike_counts

    all_pcs = np.zeros((0, pc_features.shape[1], channels_to_use.size))
    all_labels = np.zeros((0,))

    for idx2, cluster_id2 in enumerate(units_in_range):

        subsample = int(relative_counts[idx2])

        pcs = get_unit_pcs(cluster_id2, spike_clusters, spike_templates,
                           pc_feature_ind, pc_features, channels_to_use,
                           subsample)

        if pcs is not None and len(pcs.shape) == 3:

            labels = np.ones((pcs.shape[0],)) * cluster_id2

            all_pcs = np.concatenate((all_pcs, pcs),0)
            all_labels = np.concatenate((all_labels, labels),0)

    all_pcs = np.reshape(all_pcs, (all_pcs.shape[0], pc_features.shape[1]*channels_to_use.size))
    if ((all_pcs.shape[0] > 10)
            and not (all_labels == cluster_id).all()  # Not all labels are this cluster
            and (sum(all_labels == cluster_id) > 20)  # No fewer than 20 spikes in this cluster
            and (len(channels_to_use) > 0)):
        isolation_distance, l_ratio = mahalanobis_metrics(all_pcs, all_labels, cluster_id)

        d_prime = lda_metrics(all_pcs, all_labels, cluster_id)

        nn_hit_rate, nn_miss_rate = nearest_neighbors_metrics(all_pcs, all_labels,
                                                                             cluster_id,
                                                                             max_spikes_for_nn,
                                                                             n_neighbors)
    else:  # Too few spikes or cluster doesnt exist
        isolation_distance = np.nan
        d_prime = np.nan
        nn_miss_rate = np.nan
        nn_hit_rate = np.nan
        l_ratio = np.nan
    return isolation_distance, d_prime, nn_miss_rate, nn_hit_rate, l_ratio            
