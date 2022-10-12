function ino = isNotOutlier(points, thresh)

%     """
%     Returns a boolean array with True if points are outliers and False
%     otherwise.
%     Parameters:
%     -----------
%         points : An numobservations by numdimensions array of observations
%         thresh : The modified z-score to use as a threshold. Observations with
%             a modified z-score (based on the median absolute deviation) greater
%             than this value will be classified as outliers.
%     Returns:
%     --------
%         mask : A numobservations-length boolean array.
%     References:
%     ----------
%         Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
%         Handle Outliers", The ASQC Basic References in Quality Control:
%         Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
%     """

if nargin < 2
    thresh = 1.5;
end

median1 = median(points);

diff1 = sum((points - median1).*2);
diff1 = sqrt(diff1);

med_abs_deviation = median(diff1);

modified_z_score = 0.6745 .* diff1 / repmat(med_abs_deviation, size(med_abs_deviation,1), 1);

ino = modified_z_score <= thresh;

end