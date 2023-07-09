% Demo script for displaying the functionality of myStitch function.

close all;
%clc;

% Read the images and convert them to grayscale if not already
im1 = imread('im1.png');
im2 = imread('im2.png');

% Call myStitch function to create Stitched image
stitched_image = myStitch(im1, im2);

% Show stitched image
figure, clf;
imshow(stitched_image);
title('Stitched image');
