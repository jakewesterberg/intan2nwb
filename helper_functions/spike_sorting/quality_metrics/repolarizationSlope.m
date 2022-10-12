function repolarization_slope = repolarizationSlope(waveform, timestamps, window=20):
    
%     """ 
%     Spike repolarization slope (after maximum deflection point)
%     Inputs:
%     ------
%     waveform : numpy.ndarray (N samples)
%     timestamps : numpy.ndarray (N samples)
%     window : int
%         Window (in samples) for linear regression
%     Outputs:
%     --------
%     repolarization_slope : slope of return to baseline (V / s)
%     """

if nargin < 3
    window = 20;
end

[~,max_point] = max(abs(waveform));

waveform = -1.*waveform * (sign(waveform(max_point)));

repolarization_slope = timestamps(max_point:max_point+window) \ waveform(max_point:max_point+window);
repolarization_slope = repolarization_slope * 1e-6;

end