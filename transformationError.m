function error = transformationError(x, p1, p2, p3, p4)
% Auxiliary function for computeTransformation function. This funtion
% calculates the transformation error between predicted and actual points.
%
% Inputs:
%   x - A vector containing the rotation angle (theta) and distances (dx, dy)
%   p1 - The original position of the first point (x,y)
%   p2 - The transformed position of the first point (x,y)
%   p3 - The original position of the second point (x,y)
%   p4 - The transformed position of the second point (x,y)
%
% Output:
%   error - A vector containing the differences between the predicted and
%       actual transformed positions of p2 and p4

    theta = x(1);
    dx = x(2);
    dy = x(3);
    
    predicted_p2 = [cos(theta) -sin(theta); sin(theta) cos(theta)] * [p1(1); p1(2)] + [dx; dy];
    predicted_p4 = [cos(theta) -sin(theta); sin(theta) cos(theta)] * [p3(1); p3(2)] + [dx; dy];
    
    error = [predicted_p2 - p2; predicted_p4 - p4];
end