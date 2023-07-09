function [inlier_matchingPoints, outlier_matchingPoints] = findInlierOutlier(matchedPointsAll, r, H)
% Function to classify the matching points as inliers or outliers, based on
% a threshold distance (r)
%
% Inputs:
%   matchedPointsAll - a matrix containing the matching points
%   r - the threshold distance for classifying a point as an inlier or
%       outlier. Points with distance less than r are classified as inliers.
%   H - the transformation matrix.
%
% Outputs:
%   inlier_matchingPoints - a vectror containing the indexes to
%       matchedPointsAll of the inlier point pairs.
%   outlier_matchingPoints - a vectore containing the indexes to
%       matchingPointsAll of the outlier point pairs

    % Export information from arguments
    n = size(matchedPointsAll, 2);
    inlier_matchingPoints = zeros(n, 1);
    outlier_matchingPoints = zeros(n, 1);

    inliers_found = 0;
    outliers_found = 0;

    theta = H{1};
    d = H{2};

    T = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    for i = 1:n
        % Calculate norm between transformed p1 and actual p2
        p1 = matchedPointsAll(1:2, i);
        p2 = matchedPointsAll(3:4, i);
    
        expected_p2 = round(T*p1 + d);
        distance = norm(p2 - expected_p2);

        % Inlier pair of points
        if distance < r
            % find and export a pair with zero distance for stitching
            inliers_found = inliers_found + 1;
            inlier_matchingPoints(inliers_found) = i;
        else
            outliers_found = outliers_found + 1;
            outlier_matchingPoints(outliers_found) = i;
        end
    end

    % Keep only the usefull part of the two vectors
    inlier_matchingPoints = inlier_matchingPoints(1:inliers_found);
    outlier_matchingPoints = outlier_matchingPoints(1:outliers_found);

end