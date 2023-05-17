function [training_1, testing_1, training_2, testing_2, training_3, testing_3] = splitDataset(dataset_1, dataset_2, dataset_3)
%SPLITDATASET Summary of this function goes here
%   Detailed explanation goes here

    number_d1 = length(dataset_1);
    number_d2 = length(dataset_2);
    number_d3 = length(dataset_3);

    idx1 = randperm(number_d1);
    idx2 = randperm(number_d2);
    idx3 = randperm(number_d3);

    split1 = round(0.7 * number_d1);
    split2 = round(0.7 * number_d2);
    split3 = round(0.7 * number_d3);

    training_1 = dataset_1(idx1(1:split1), :);
    testing_1 = dataset_1(idx1(split1+1:end), :);

    training_2 = dataset_2(idx2(1:split2), :);
    testing_2 = dataset_2(idx2(split2+1:end), :) ;

    training_3 = dataset_3(idx3(1:split3), :);
    testing_3 = dataset_3(idx3(split3+1:end), :);

end