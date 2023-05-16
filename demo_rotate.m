% Demo script for demostrating the use of rotateImage function.

close all;
clear;

% Read the image
img = imread('images/text1.png');

% Convert the image to grayscale 
gray_img = rgb2gray(img);

% Create a figure and display the grayscale image to subplot 1
figure('Name', 'demo_rotate'), clf, hold on;
subplot(1,2,1);
imshow(gray_img);
title('Input grayscale image')

% Rotate the image, using the rotateImage function
angle = 60; %deg
rotated_img = rotateImage(gray_img, angle);

% Save the rotated image
imwrite(rotated_img, 'images/rotated_60.png')
% Display the rotated image
subplot(1,2,2)
imshow(rotated_img);
title('Rotated image');