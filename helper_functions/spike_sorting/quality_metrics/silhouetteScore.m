function ss = silhouetteScore(spike_clusters,spike_templates,total_units,pc_features,pc_feature_ind,total_spikes)

%     def score_inner_loop(i, cluster_ids):
%         """
%         Helper to loop over cluster_ids in one dimension. We dont want to loop over both dimensions in parallel-
%         that will create too much worker overhead
%         Args:
%             i: index of first dimension
%             cluster_ids: iterable of cluster ids
%         Returns: scores for dimension j
%         """

cluster_ids = unique(spike_clusters);

random_spike_inds = randperm(numel(spike_clusters));
random_spike_inds = random_spike_inds(1:total_spikes);
num_pc_features = size(pc_features,2);
num_channels = max(pc_feature_ind) + 1;

all_pcs = zeros(total_spikes, num_channels * num_pc_features);

idx = 1;
for i = random_spike_inds

    unit_id = spike_templates(i);
    channels = pc_feature_ind(unit_id,:);

    for j = 1:num_pc_features
        all_pcs(idx, channels + num_channels * j) = pc_features(i,j,:);
    end
    idx = idx + 1;
end

cluster_labels = squeeze(spike_clusters(random_spike_inds));

SS = nan(total_units, total_units);

scores = [score_inner_loop(i, cluster_ids) 
    for i = cluster_ids

    for i, col_score in enumerate(scores):
        for j, one_score in enumerate(col_score):
            SS[i, j] = one_score
        end
    end

    a = np.nanmin(SS, 0)
    b = np.nanmin(SS, 1)

    ss=  np.array([np.nanmin([a,b]) for a, b in zip(a,b)])

    end


function score_inner_loop(i, cluster_ids)
scores_1d = [];
for j = cluster_ids
    if j > i
        inds = np.in1d(cluster_labels, np.array([i, j]))
        X = all_pcs[inds, :]
        labels = cluster_labels[inds]

        # len(np.unique(labels))=1 Can happen if total_spikes is low:
        if (len(labels) > 2) and (len(np.unique(labe gls)) > 1):
            scores_1d.append(silhouette_score(X, labels))
        else:
            scores_1d.append(np.nan)
        end
    else:
        scores_1d.append(np.nan)
    end
end
