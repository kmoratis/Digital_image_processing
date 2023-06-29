function [rotatedImage, rotatedP] = rotateImage(I, angle, p, show_im)
%ROTATEIMAGE Summary of this function goes here
%   Detailed explanation goes here

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