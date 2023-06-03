% Demo script for demonstrating the use of createModel and readText
% functions
close all;

% Create k-NN models
N = [5 5 5];
seed = 3;
K = [5 1 1];
[knnModel_1, knnModel_2, knnModel_3] = createModel('text1_v3.png', 'text1_v3.txt', N, seed, K);

% Call readText function
%img_path = 'my.png';
img_path = 'text2_150dpi_rot.png';
txt_path = 'text2.txt';
readText(img_path, txt_path, knnModel_1, knnModel_2, knnModel_3, N);
