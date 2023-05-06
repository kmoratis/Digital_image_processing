close all;
clear;
clc;

% read image
img = imread('images/rotated_45.png');
%img = imread('images/text1.png');
%img = rgb2gray(img);

% find first angle evaluation, using 2D DFT
angle = findRotationAngle(img);