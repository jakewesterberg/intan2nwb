function out_mat = pullVecs(in_mat, idxs, vrange)

nd = ndims(in_mat);
out_mat = [];
for ii = 1 : numel(idxs)
    out_mat = cat(nd+1, out_mat, in_mat(:,idxs(ii)-vrange(1):idxs(ii)+vrange(2)));
end

end