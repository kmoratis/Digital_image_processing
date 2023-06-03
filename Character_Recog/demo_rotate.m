% Demo script for demostrating the use of rotateImage function.

close all;
clear;

% Read the image
img = imread('text1_v3.png');

% Rotation angle
angle = 60; %deg

% Convert the image to grayscale 
gray_img = rgb2gray(img);

[rows, cols] = size(gray_img);

% Add one pixel white-padding to each side
temp = uint8(zeros(rows+2, cols+2) + 255);
temp(2:rows+1, 2:cols+1) = gray_img(:,:);
gray_img = temp;

% Create a figure and display the grayscale image to subplot 1
figure('Name', 'demo_rotate'), clf, hold on;
subplot(1,2,1);
imshow(gray_img);
title('Input grayscale image')

% Rotate the image, using the rotateImage function
rotated_img = rotateImage(gray_img, angle);

% Save the rotated image
%imwrite(rotated_img, 'rotated_V3_60.png')

% Display the rotated image
subplot(1,2,2)
imshow(rotated_img);
title('Rotated image');
