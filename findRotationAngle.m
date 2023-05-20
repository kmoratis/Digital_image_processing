function angle = findRotationAngle(x, precision, fig_show)
% Finds a more precise evaluation of the rotation angle of the image, by
% doing a linear search in an area near the evaluated angle given by
% angleEvaluation function.
%
% Inputs:
%   x - the grayscale image, which rotation angle we want to measure
%   presicion - number of points to be examined in linear search
%   fig_show - boolean variable, indicating whether the figures of
%   angleEvaluation will be shown
%
% Output:
%   angle - an evaluation of the rotation angle of the input image

    % Call the angleEvaluation function to find a first evaluation of the
    % angle
    angle_eval = angleEvaluation(x, fig_show);

    linear_num_points = precision;

    % Linear search in a region near the evaluated value
    candidate_angles = linspace(angle_eval-1, angle_eval+1, linear_num_points);

    candidate_angles = [angle_eval candidate_angles];

    % Initialize the variables
    maxVar = -Inf;
    angle = 0;

    % Loop through each candidate angle
    for i = 1:length(candidate_angles)
        % Undo the rotation, using function rotateImage
        unrotated_img = rotateImage(x, -candidate_angles(i));

        % Calculate the brightness projection to the vertical axis
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
