function wavePTR(waveform):

    """ 
    Peak-to-trough ratio of 1D waveform
    Inputs:
    ------
    waveform : numpy.ndarray (N samples)
    Outputs:
    --------
    PT_ratio : waveform peak-to-trough ratio
    """

    trough_idx = np.argmin(waveform)

    peak_idx = np.argmax(waveform)

    PT_ratio = np.abs(waveform[peak_idx] / waveform[trough_idx])

    return PT_ratio
