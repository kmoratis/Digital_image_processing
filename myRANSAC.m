function [H, inlierMatchingPoints, outlierMatchingPoints] = myRANSAC(matchingPoints, r, N)
% Function implementing the RANSAC algorithm for estimating a
% transformation between two sets of matching points.
%
% Inputs:
%   matchingPoints - a matrix containing the matching points. Each column
%       represents a pair of points. The first two rows contain the (x,y)
%       of the first pixel and the second two contain the (x,y) of its
%       match.
%   r - the threshold distance for inlier/ outlier classification
%   N - the number of iterations 
%
% Outputs:
%   H - a cell array, containing the {theta, d} of the transformation
%   inlierMatchingPoints - a vector containing the inlier matching points.
%       Each value contains an index to a matchingPoints pair. 
%   outlierMatchingPoints - a vector containing the outlier matching
%       points. Each value contains an index to a matchingPoints pair.
    
    % RANSAC algorithm implementation
    max_score = 1e+10;
    best_H = 0;
    n = size(matchingPoints, 2);
    
    % Create a vector containing a random permutation of the integers 1-n
    p = randperm(n);
    
    for i=1:N
        % Pick two random pairs.
        idx1 = p(i);
        idx2 = p(i+1);
        sample(:,1) = matchingPoints(:, idx1);
        sample(:,2) = matchingPoints(:, idx2);
    
        % Compute the transformation, defined by the two pairs selected
        H = computeTransformation(sample);
        
        % Calculate the score of the candidate H
        score = computeScore(matchingPoints, H);
    
        % Keep the min score ( min Euclidean distance sum ).
        if score < max_score
            max_score = score;
            best_H = H;
        end
    end
    
    % With the best_H transform and the specific r, we will find the
    % inlier and outlier matching points
    [inlierMatchingPoints, outlierMatchingPoints] = findInlierOutlier(matchingPoints, r, best_H);

    H = best_H; % return best_H
end