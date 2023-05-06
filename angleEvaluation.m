function angle = angleEvaluation(img)
% Gives a first evaluation of the rotation angle of the input image, using
% the 2D DFT. The angle is measured from the vertical axis, in a counter-clockwise direction.
% Notes:
%   - the function works for angles in the range [0, 180), due to the
% symmetry of the DFT.
%   - the function is not exactly accurate, because of the limited analysis
%   of the DFT
%
% Input:
%   img - the image in grayscale which rotation angle we want to measure
%
% Output:
%   angle - the evaluation of the angle

    % Blur the image using Gaussian kernel
    blur_kernel = fspecial('gaussian', [10 15], 15);
    blurredImage = conv2(img, blur_kernel, 'same');
    
    figure('Name', 'Angle Evaluation');
    subplot(2,2,[1,2]);
    imshow(blurredImage, []);
    title('Blurred image');

    % Compute the 2D DFT of the image
    dft = fft2(blurredImage);

    % Remove dc value
    dft(1,1) = 0;
    
    % Shift the zero-frequency component of the DFT to the center
    dft = fftshift(dft);

    % Compute the magnitude of the DFT
    mag_spectrum = abs(dft);
    mag_spectrum = 10 * log(1 + mag_spectrum);
    mag_spectrum = max(mag_spectrum, 0);
    subplot(2,2,3);
    imshow(mag_spectrum, []);
    title('Magnitude spectrum');
    
    % Find the peak of the magnitude of the DFT
    [M, ~] = max(mag_spectrum(:));

    % Create a threshold value
    threshold = 0.9*M;

    % Binarize image, using the threshold
    bin_img = imbinarize(mag_spectrum, threshold);
    subplot(2,2,4);
    imshow(bin_img);
    title('Binarized image');

    % Find points with magnitude more that the treshold
    [row, col, ~] = find(bin_img);

    % Find the index of the most left and upper pixel
    point_idxs = row + col;
    [~, idx] = min(point_idxs);

    % Find image center coordinates
    center_x = size(mag_spectrum, 2)/2 - 1;
    center_y = size(mag_spectrum, 1)/2 - 1;

    % Find the coordinates of the most left non-zero point in the binarized picture 
    y = row(idx);
    x = col(idx);

    % Find the cos of the angle between upper left point and center point
    dist = sqrt((y-center_y)^2 + (center_x-x)^2);
    cos_angle = (center_y-y) / dist;

    angle = acosd(cos_angle);

    % Display the result
    disp(['First evaluation using DFT: The image is rotated by ', num2str(angle), ' degrees.']);

end
