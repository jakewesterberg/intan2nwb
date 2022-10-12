function PT_ratio = wavePTR(waveform)

%     """ 
%     Peak-to-trough ratio of 1D waveform
%     Inputs:
%     ------
%     waveform : numpy.ndarray (N samples)
%     Outputs:
%     --------
%     PT_ratio : waveform peak-to-trough ratio
%     """

    [~,trough_idx] = min(waveform);
    [~,peak_idx] = max(waveform);
    PT_ratio = abs(waveform(peak_idx) / waveform(trough_idx));

end
