% Test for matching descriptors

close all;
clc;

% Read the images and convert them to grayscale if not already
im1 = imread('im1.png');

im2 = imread('im2.png');

% Get the size of the image
[~, ~, numberOfColors] = size(im1);

% If the image is RGB, convert it to grayscale
if numberOfColors > 1
    im1 = rgb2gray(im1);
end  

% Get the size of the image
[~, ~, numberOfColors] = size(im2);

% If the image is RGB, convert it to grayscale
if numberOfColors > 1
    im2 = rgb2gray(im2);
end

% Harris Detector parameters
sigma = 0.33; % Gaussian smoothing parameter %0.3
threshold = 5000; % Threshold
k = 0.22; % Harris k constant %0.18
show_im = false; % Parameter for img display

% Detect corner points for image 1
points_1 = myDetectHarrisFeatures(im1, sigma, k, threshold, show_im);

% Detect corner points for image 2
points_2 = myDetectHarrisFeatures(im2, sigma, k, threshold, show_im);

% Perform the descriptor matching, by calling the appropriate function
percentageThreshold = 0.8; % threshold for keeping only percentage of the points with the strongest simillarity
descriptor_upgrade = false; % which version of the descriptor to use
show_im_matching = true; % whether to show the images with the matching points or not
descriptorMatching(im1, im2, points_1, points_2, percentageThreshold, descriptor_upgrade, show_im_matching);