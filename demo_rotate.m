close all;
clear;

% read image
img = imread('text1.png');

% convert image to grayscale 
gray_img = rgb2gray(img);

% find angle evaluation
%angle_eval = angleEvaluation(gray_img);

% rotate image
angle = 20; %deg
interp_method = 'bilinear'; 
rotated_img = rotateImage(gray_img, angle, interp_method);
% save rotated image
imwrite(rotated_img, 'rotated_20.png')
% display rotated image
figure();
imshow(rotated_img);
title('Rotated Image');