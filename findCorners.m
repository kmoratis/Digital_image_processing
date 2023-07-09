function c = findCorners(I, sigma, k, Rthres)
% Function finding the corner points of an input image, using the Harris
% corner detection algorithm.
%
% Inputs:
%   I - the grayscale image
%   sigma - a float parameter, used for Gaussian mask. 
%       ( Bigger sigma values lead to slower reduction of the Gaussian
%       mask, which leads to more neighbor pixels taken into account. )
%   k - a float parameter for Harris detection algorithm. 
%       ( Bigger k values, makes the algorithm more biased towards
%       detecting corners in contrast to edges. )
%   Rthres - threshold value for keeping only the features (corners, edges)
%       with the strongest R value. (most likely to be true positive)
%
% Output:
%   c - A boolean vector, containing 1 if the particular pixel of the image
%       is a corner point.

    % Derivative masks
    [dx, dy] = meshgrid(-1:1, -1:1);

    % Compute the image derivatives
    Ix = conv2(double(I), dx, 'same');
    Iy = conv2(double(I), dy, 'same');

    % Compute products
    Ix2 = Ix.^2;
    Iy2 = Iy.^2;
    Ixy = Ix.*Iy;

    % Compute the sums of products at each pixel
    Sx2 = conv2(Ix2, fspecial('gaussian', 2*ceil(3*sigma)+1, sigma), 'same');
    Sy2 = conv2(Iy2, fspecial('gaussian', 2*ceil(3*sigma)+1, sigma), 'same');
    Sxy = conv2(Ixy, fspecial('gaussian', 2*ceil(3*sigma)+1, sigma), 'same');

    % Compute the Harris response
    R = (Sx2.*Sy2 - Sxy.^2) - k*(Sx2 + Sy2).^2;

    % Check if pixel is a corner
    c = R > Rthres;
end