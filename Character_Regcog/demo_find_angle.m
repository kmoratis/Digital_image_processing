% Demo script for demostrating the use of findRotationAngle function.

close all;
clear;
clc;

path = "text1_rot.png";

% Read the image
img = imread(path);

% Get the size of the image
[rows, cols, numberOfColors] = size(img);

% If the image is RGB, convert it to grayscale
if numberOfColors > 1
    img = rgb2gray(img);
end

% Preprocessing to remove the black padding
if strcmp(path, "text2_150dpi_unrot.png") || strcmp(path, "text1_rot.png")
    [m, n] = size(img);
    img(1:9, :) = 255; % top
    img(m-9:m, :) = 255; % bottom
    img(:, 1:7) = 255; % left 
    img(:, n-7:n) = 255; % right
end

% Only for text1_v3
if strcmp(path, "text1_v3.png")
    % Add one pixel white-padding to each side
    temp = uint8(zeros(rows+2, cols+2) + 255);
    temp(2:rows+1, 2:cols+1) = img(:,:);
    img = temp;
    imshow(img);
end

% Find the rotation angle of the input image
precision = 15; % number of points to be examined in linear search
fig_show = true;
angle = findRotationAngle(img, precision, fig_show);

% Call rotatImage function to reverse the image rotation
unrot_img = rotateImage(img, -angle);
figure(), clf;
imshow(unrot_img);
title('Unrotated image');
