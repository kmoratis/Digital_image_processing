% Demo script for demostrating the use of splitWords and splitCharacters functions.

%clear;
%clc;
%close all;

path = "text1_v3.png";

img = imread(path);

% If in RGB, transform to grayscale
[m, n, numColors] = size(img);
if numColors==3
    gray_img = rgb2gray(img);
else
    gray_img = img;
end

% Preprocessing 
% If image is text1, add a 10px white-padding in each side
if strcmp(path, "text1_v3.png")
    temp = uint8(zeros(m+20, n+20) + 255);
    temp(11:m+10, 11:n+10) = gray_img;
    gray_img = temp;
end
% If image is text2, make edges (2px each) white
if strcmp(path, "text2_150dpi_rot.png")
    
    % Preprocessing: If image is text2, make edges (2px each) white
    gray_img(1:2, :) = 255;
    gray_img(m-1:m, :) = 255;
    gray_img(:, 1:2) = 255;
    gray_img(:, n-1:n) = 255;
end

% Call splitWords function
words = splitWords(gray_img);
disp(['Words detected: ', num2str(length(words))]);

% Call splitCharacters function
chars = splitCharacters(words);
disp(['Characters detected: ', num2str(length(chars))]);
