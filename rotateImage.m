function [rotatedImage, rotatedP] = rotateImage(I, angle, p, show_im)
% Auxiliary function for image rotation, and for finding pixel's coordinates in
% the rotated image.
%
% Inputs:
%   I - a grayscale image, which we want to rotate
%   angle - a number specifying the rotation angle (in degrees)
%   p - a pixel (x,y) in the original image (I), which coordinates we want
%       to find in the transformed one
%   show_im - a boolean parameter, specifying whether the function will
%       show some relative to the calculations images
%
% Ouputs:
%   rotatedImage - the grayscale input image (I), rotated by (angle) degrees.
%   rotatedP - the coordinates of the input pixel (p) in the rotatedImage.

    p = p';

    % Image rotation
    rotatedImage = imrotate(I, angle, 'crop');

    rotMatrix = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];

    % Calculate the center of the input image
    imCenter = (size(I, 1, 2)/2)';
    % Calc the center of the rotater image
    rotImCenter = (size(rotatedImage, 1, 2)/2)';

    rotatedP = rotMatrix * (p - imCenter) + rotImCenter;

    rotatedP = rotatedP';

    % If show_im == true, show the two images with points
    if show_im
        figure, clf;
        imshow(I);
        hold on
        plot(p(2), p(1), 'r+', 'LineWidth', 20);
        plot(imCenter(2), imCenter(1), 'g+', 'LineWidth', 20);
        title('Original Image');

        figure, clf;
        imshow(rotatedImage);
        hold on
        plot(rotatedP(2), rotatedP(1), 'r+', 'LineWidth', 20);
        plot(rotImCenter(2), rotImCenter(1), 'g+', 'LineWidth', 20);
        title('Rotated Image');
    end

end