function fraction_missing = amplitudeCutoff(amplitudes, num_histogram_bins = 500, histogram_smoothing_value = 3)

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

h = histcounts(amplitudes, num_histogram_bins, density=True);

pdf = gaussian_filter1d(h,histogram_smoothing_value)
support = b[:-1]

peak_index = np.argmax(pdf)
G = np.argmin(np.abs(pdf[peak_index:] - pdf[0])) + peak_index

bin_size = np.mean(np.diff(support))

fraction_missing = sum(pdf[G:])*bin_size
fraction_missing = min([fraction_missing, 0.5]);

end