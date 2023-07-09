function H = computeTransformation(sample)
% Function calculating the transformation, given two pairs of pixels, using
% the MSE estimation.
%
% Inputs: 
%   sample - 4x2 array containing [pair1 pair2] vectors, 
%       where pair1 = [x1 y1 x2 y2]' and pair2 = [x3 y3 x4 y4]'.
%
% Output:
%   H - 1x2 cell array, containing {theta, d} of the transformation.

    % Extract the coordinates of the points
    % p1, p2 pair
    p1 = sample(1:2, 1);
    p2 = sample(3:4, 1);
    
    % p3, p4 pair
    p3 = sample(1:2, 2);
    p4 = sample(3:4, 2);

    % Initial estimations for theta, d1, d2
    x0 = [0; 0; 0];

    % Compute the rotation angle (theta) and translation vector (d), using
    % the lsqnonlin function -> MSE estimation

    x = lsqnonlin(@(x) transformationError(x, p1, p2, p3, p4), x0);

    % Retrieve values
    theta = x(1);
    d = x(2:3);

    H = {theta, d};
end