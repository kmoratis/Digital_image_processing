% Demo script for demostrating the use of findRotationAngle function.

close all;
clear;
clc;

% Read the image
img = imread('images/rotated_45.png');

% Get the size of the image
[rows, cols, numberOfColors] = size(img);

% If the image is RGB, convert it to grayscale
if numberOfColors > 1
    img = rgb2gray(img);
end

% Find the rotation angle of the input image
precision = 15; % number of points to be examined in linear search
fig_show = true;
angle = findRotationAngle(img, precision, fig_show);