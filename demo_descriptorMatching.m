% Demo script for demonstrating the use of descriptorMatching function.

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

% Descriptor parameters
rhom = 2; % Min rad
rhoM = 12; % Max rad
rhostep = 1; % Rad step
N = 40; % Number of points

desc_upgrade = false; % which version of the descriptor to use

% Calc descriptor vector size
s = 1 + (rhoM - rhom) / rhostep;

% Create array to store actual corner points ( descriptor not empty )
cor_1 = zeros(size(points_1));
cor1_found = 0;

% Create array to store im1 descriptors
N1 = size(points_1, 1);
desc1 = zeros(s, N1);

% Calculate im1 corners descriptors
for i = 1:N1
    p = points_1(i, :);

    if desc_upgrade == true % use desc_upgrade verion of the descriptor
        d = myLocalDescriptorUpgrade(im1, p, rhom, rhoM, rhostep, N);
        if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
            cor1_found = cor1_found + 1;
            desc1(:, cor1_found) = d(:);
            cor_1(cor1_found, :) = p;
        end
    else
        d = myLocalDescriptor(im1, p, rhom, rhoM, rhostep, N);
        if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
            cor1_found = cor1_found + 1;
            desc1(:, cor1_found) = d(:);
            cor_1(cor1_found, :) = p;
        end
    end
end

% keep the usefull part of the array
cor_1 = cor_1(1:cor1_found, :);
desc1 = desc1(:, 1:cor1_found);

% Create array to store actual corner points ( descriptor not empty )
cor_2 = zeros(size(points_2));
cor2_found = 0;

% Create array to store im2 descriptors
N2 = size(points_2, 1);
desc2 = zeros(s, N2);

% Calculate im2 corners descriptors
for i = 1:N2
    p = points_2(i, :);

    if desc_upgrade == true % use desc_upgrade verion of the descriptor
        d = myLocalDescriptorUpgrade(im2, p, rhom, rhoM, rhostep, N);
        if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
            cor2_found = cor2_found + 1;
            desc2(:, cor2_found) = d(:);
            cor_2(cor2_found, :) = p(1, :);
        end
    else
        d = myLocalDescriptor(im2, p, rhom, rhoM, rhostep, N);
        if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
            cor2_found = cor2_found + 1;
            desc2(:, cor2_found) = d(:);
            cor_2(cor2_found, :) = p;
        end
    end
end

% keep the usefull part of the array
cor_2 = cor_2(1:cor2_found, :);
desc2 = desc2(:, 1:cor2_found);

% Perform the descriptor matching, by calling the appropriate function
percentageThreshold = 0.8; % threshold for keeping only percentage of the points with the strongest simillarity
show_im_matching = true; % whether to show the images with the matching points or not

matchingPoints = descriptorMatching(desc1, desc2, percentageThreshold);

if show_im_matching
    colors = ['r', 'y', 'b', 'c', 'm'];
    figure, clf;
    % Image 1
    subplot(1,2,1);
    imshow(im1);
    hold on
    for i = 1:5
        % Find best {i} matches and plot them 
        p1_idx = matchingPoints(1,i); % best_i -> index of best point of img1
        pixel1 = cor_1(p1_idx, :); % [x, y] in image 1
        
        % Display the corners on the original image
        plot(pixel1(1, 2), pixel1(1, 1), append(colors(i), '+'), 'LineWidth', 10);
        
    end
    hold off
    % Image 2
    subplot(1,2,2);
    imshow(im2);
    hold on
    for i = 1:5
        p2_idx = matchingPoints(2,i); % best_j
        pixel2 = cor_2(p2_idx, :); % [x, y] in image 2
        
        plot(pixel2(1,2), pixel2(1,1), append(colors(i), '+'), 'LineWidth', 10);
        
    end
    hold off
end