function [velocity_above, velocity_below] = getVelocity(channels, times, distance_between_channels)

%     """
%     Calculate slope of trough time above and below soma.
%     Inputs:
%     -------
%     channels : np.ndarray
%         Channel index relative to soma
%     times : np.ndarray
%         Trough time relative to peak channel
%     distance_between_channels : float
%         Distance between channels (m)
%     Outputs:
%     --------
%     velocity_above : float
%         Inverse of velocity of spike propagation above the soma (s / m)
%     velocity_below : float
%         Inverse of velocity of spike propagation below the soma (s / m)
%     """

if nargin < 3
    distance_between_channels = 10e-6;
end

above_soma = channels >= 0;
below_soma = channels <= 0;

if np.sum(above_soma) > 1
    slope_above = channels(above_soma) \ times(above_soma);
    velocity_above = slope_above / distance_between_channels;
else
    velocity_above = NaN;
end

if np.sum(below_soma) > 1
    slope_below = channels(below_soma) \ times(below_soma);
    velocity_below = slope_below / distance_between_channels;
else
    velocity_below = NaN;
end

end