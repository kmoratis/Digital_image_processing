function angle = findRotationAngle(x)
% Finds a more precise evaluation of the rotation angle of the image, by
% doing a linear search in an area near the evaluated angle given by
% angleEvaluation funtion.
%
% Input:
%   x - the grayscale image, which rotation angle we want to measure
%
% Output:
%   angle - an evaluation of the rotation angle of the input image

    % call angleEvaluation function to find a first evaluation of the angle
    angle_eval = angleEvaluation(x);

    % linear search in a region near evaluated value
    candidate_angles = linspace(angle_eval-5, angle_eval+5, 100);

    % Initialize variables
    maxVar = -Inf;
    angle = 0;
    interp_method = 'bilinear';

    % Loop through each candidate angle
    for i = 1:length(candidate_angles)
        % Undo the rotation, using function rotateImage
        unrotated_img = rotateImage(x, -candidate_angles(i), interp_method);

        % Calculate brightness projection to the vertical axis
        brightness = sum(unrotated_img, 2);

        % Compute the variation in brightness
        var = max(brightness) - min(brightness);

        % Check if the current variation is bigger than the max found
        if var > maxVar
            maxVar = var;
            angle = candidate_angles(i);
        end
    end

    % Display the result
    disp(['After linear search: The image is rotated by ', num2str(angle), ' degrees.']);

end
