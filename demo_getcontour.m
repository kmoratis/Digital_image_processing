% Demo script for demostrating the use of getcontour function.

clear;
clc;
close all;

% Read the image containing the letter
img = imread('Letters/letter_e.png');

% Convert to grayscale
gray_img = rgb2gray(img);

% Call getcontour function
fig_shown = true;
c = getcontour(gray_img, fig_shown);

% Call display contour function to display the contours found
[m, n] = size(gray_img);
displaycontour(c, m ,n);