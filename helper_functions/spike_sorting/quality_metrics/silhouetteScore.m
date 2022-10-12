function silhouetteScore(spike_clusters,spike_templates,total_units,pc_features,pc_feature_ind,total_spikes):

%     def score_inner_loop(i, cluster_ids):
%         """
%         Helper to loop over cluster_ids in one dimension. We dont want to loop over both dimensions in parallel-
%         that will create too much worker overhead
%         Args:
%             i: index of first dimension
%             cluster_ids: iterable of cluster ids
%         Returns: scores for dimension j
%         """

        scores_1d = []
        for j in cluster_ids:
            if j > i:
                inds = np.in1d(cluster_labels, np.array([i, j]))
                X = all_pcs[inds, :]
                labels = cluster_labels[inds]

                # len(np.unique(labels))=1 Can happen if total_spikes is low:
                if (len(labels) > 2) and (len(np.unique(labels)) > 1):
                    scores_1d.append(silhouette_score(X, labels))
                else:
                    scores_1d.append(np.nan)
            else:
                scores_1d.append(np.nan)
        return scores_1d

    cluster_ids = np.unique(spike_clusters)

    random_spike_inds = np.random.permutation(spike_clusters.size)
    random_spike_inds = random_spike_inds[:total_spikes]
    num_pc_features = pc_features.shape[1]
    num_channels = np.max(pc_feature_ind) + 1

    all_pcs = np.zeros((total_spikes, num_channels * num_pc_features))

    for idx, i in enumerate(random_spike_inds):

        unit_id = spike_templates[i]
        channels = pc_feature_ind[unit_id,:]

        for j in range(0,num_pc_features):
            all_pcs[idx, channels + num_channels * j] = pc_features[i,j,:]

    cluster_labels = np.squeeze(spike_clusters[random_spike_inds])

    SS = np.empty((total_units, total_units))
    SS[:] = np.nan


    if do_parallel:
        from joblib import Parallel, delayed
        scores = Parallel(n_jobs=-1, verbose=2)(delayed(score_inner_loop)(i, cluster_ids) for i in cluster_ids)
    else:
        scores = [score_inner_loop(i, cluster_ids) for i in cluster_ids]

    # Fill the 2d array
    for i, col_score in enumerate(scores):
        for j, one_score in enumerate(col_score):
            SS[i, j] = one_score

    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      a = np.nanmin(SS, 0)
      b = np.nanmin(SS, 1)

    return np.array([np.nanmin([a,b]) for a, b in zip(a,b)])