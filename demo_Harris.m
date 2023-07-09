% Demo script for demonstrating the use of myDetectHarrisFeatures and findCorners functions.

close all;
clc;

% Read the image
I = imread('im1.png');
%I = imread('imForest1.png');

% Get the size of the image
[height, width, numberOfColors] = size(I);

% If the image is RGB, convert it to grayscale
if numberOfColors > 1
    I = rgb2gray(I);
end

% Gaussian smoothing parameter
sigma = 0.9; %2.2

% Corner threshold
threshold = 5e+6; 

% Harris k constant
k = 0.23;

% Parameter for img display
show_im  = true;

% Detect corner points, using Harris Detection Algorithm and display them
corners = myDetectHarrisFeatures(I, sigma, k, threshold, show_im);
