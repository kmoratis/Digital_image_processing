function y = rotateImage(x, angle, interp_method)
% Rotate the input image x by the specified angle (in degrees) counter-clockwise, with white
% padding 
%
% Inputs:
%   x - the input image to be rotated
%   angle - the angle of rotation in degrees
%   interp_method - the interpolation method for the imrotate
%
% Output:
%   y - the rotated image with white padding

    % compute negative of image x
    neg_x = 255 - x;

    % rotate the negative x image
    neg_y = imrotate(neg_x, angle, interp_method);

    % compute the negative of image neg_y
    y = 255 - neg_y;
end
