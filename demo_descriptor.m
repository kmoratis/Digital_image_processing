% Demo script for demonstrating the use of myLocalDescriptor  and 
% myLocalDescriptorUpgrade functions.

close all;
%clc;

% Read the image
I = imread('im1.png');

% Get the size of the image
[~, ~, numberOfColors] = size(I);

% If the image is RGB, convert it to grayscale
if numberOfColors > 1
    I = rgb2gray(I);
end

rhom = 5; %5
rhoM = 10; %10
rhostep = 1;
N = 40; %40

p = [500 500];

% Find the descriptor of the p
d = myLocalDescriptor(I, p, rhom, rhoM, rhostep, N);

% Find rotated image and new coordinates of p
angle = 45;
show_im = true;
[rotatedImage, rotatedP] = rotateImage(I, angle, p, show_im);

rotatedP = round(rotatedP);

% Find the descriptor of the p in the rotated image
d_rot = myLocalDescriptor(rotatedImage, rotatedP, rhom, rhoM, rhostep, N);
