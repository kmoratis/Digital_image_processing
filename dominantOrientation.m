function dominant_orientation = dominantOrientation(Gdir_quantized, x, y)
%DOMINANTORIENTATION Summary of this function goes here
%   Detailed explanation goes here

    % Find patch around keypoint
    patch_size = 16;
    half_patch_size = patch_size / 2;
    patch = Gdir_quantized(y - half_patch_size + 1 : y + half_patch_size, ...
                           x - half_patch_size + 1 : x + half_patch_size);

    % Compute the histogram of the patch
    num_bins = 36;
    histogram = histcounts(patch(:), 0:num_bins);

    % Get the max bin value
    [~, max_bin] = max(histogram);

    % Convert the bin index to angle
    dominant_orientation = (max_bin - 1) * (360 / num_bins);
end