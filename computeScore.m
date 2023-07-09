function score = computeScore(matchedPoints, candidateH)
% Function to compute the score of a candidate transformation based on the
% difference between the expected and the actual transformed positions of
% matched points.
%
% Inputs:
%   matchedPoints - matrix containing the matched points. Each column
%       represents a pair of points, where the first two rows represents
%       the first point (x; y) and the second two represents its match.
%   candidateH - a cell array containing the candidate transformation
%       parameters {theta, d}
%
% Output:
%   score - the score of the candidate transformation.

    % Export information from arguments
    n = size(matchedPoints, 2);
    theta = candidateH{1};
    d = candidateH{2};
    T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    
    score = 0;
    
    for i = 1:n
        p1 = matchedPoints(1:2, i);
        p2 = matchedPoints(3:4, i);
    
        expected_p2 = T*p1 + d;
    
        score = score + norm(p2 - expected_p2); 
    end
end