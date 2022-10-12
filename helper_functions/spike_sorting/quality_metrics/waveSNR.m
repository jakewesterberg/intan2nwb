function wsnr = waveSNR(W)
%
%     """
%     Calculate SNR of spike waveforms.
%     Converted from Matlab by Xiaoxuan Jia
%     ref: (Nordhausen et al., 1996; Suner et al., 2005)
%     Input:
%     -------
%     W : array of N waveforms (N x samples)
%     Output:
%     snr : signal-to-noise ratio for unit (scalar)
%     """

W_bar = nanmean(W);
A = max(W_bar) - min(W_bar);
e = W - repmat(W_bar, size(W,1), 1);
wsnr = A/(2.*nanstd(e));

end