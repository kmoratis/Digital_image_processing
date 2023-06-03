% Demo script for demostrating the use creatModel functions
close all;

% Initialize variables
path = "text1_v3.png";
txt_path = 'text1_v3.txt';

N = [500 500 500]; % N1, N2, N3 values for get_contour interpolation
seed = 3;
K = [3, 1, 1]; % k1, k2, k3 values for the 3 k-NN models

% Call createModel function
createModel(path, txt_path, N, seed, K);
