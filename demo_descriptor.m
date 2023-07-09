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

rhom = 5;
rhoM = 20; %10
rhostep = 1;
N = 8; %40

p = [202 202];

% Find the descriptor of the p
d = myLocalDescriptorUpgrade(I, p, rhom, rhoM, rhostep, N); % use myLocalDescriptorUpgrade || myLocalDescriptor here
disp('Descriptor of the pixel in the given image: ');

% Find rotated image and new coordinates of p
angle = 45;
show_im = true;
[rotatedImage, rotatedP] = rotateImage(I, angle, p, show_im);

rotatedP = round(rotatedP);
disp(d');

% Find the descriptor of the p in the rotated image
d_rot = myLocalDescriptorUpgrade(rotatedImage, rotatedP, rhom, rhoM, rhostep, N); % use myLocalDescriptorUpgrade || myLocalDescriptor here
disp('Descriptor of the pixel in the rotated image: ');
disp(d_rot');
