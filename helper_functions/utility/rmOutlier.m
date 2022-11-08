function data_mat = rmOutlier(data_mat)

bad_trials = squeeze(max(abs(data_mat),[],2))';

bad_trials = squeeze(isoutlier(bad_trials, 'gesd'));

for ii = 1 : size(bad_trials, 2)
    data_mat(ii,:,bad_trials(:,ii)) = NaN;
end

end