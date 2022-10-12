function fraction_missing = amplitudeCutoff(amplitudes, num_histogram_bins, histogram_smoothing_value)

%     """ Calculate approximate fraction of spikes missing from a distribution of amplitudes
%     Assumes the amplitude histogram is symmetric (not valid in the presence of drift)
%     Inspired by metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
%     Input:
%     ------
%     amplitudes : numpy.ndarray
%         Array of amplitudes (don't need to be in physical units)
%     Output:
%     -------
%     fraction_missing : float
%         Fraction of missing spikes (0-0.5)
%         If more than 50% of spikes are missing, an accurate estimate isn't possible
%     """

if nargin < 3
    histogram_smoothing_value = 3;
end
if nargin < 4
    num_histogram_bins = 500;
end

[h, b] = histcounts(amplitudes, num_histogram_bins, 'Normalization', 'Probability');

pdf1 = gaussFilter1(1:numel(h),h,histogram_smoothing_value);
support = b(1:end-1);

[~, peak_index] = max(pdf1);
[~, G] = min(abs(pdf1(peak_index:end) - pdf1(1)));
G = G + peak_index;

bin_size = mean(diff(support));

fraction_missing = sum(pdf(G:end))*bin_size;
fraction_missing = min([fraction_missing, 0.5]);

end