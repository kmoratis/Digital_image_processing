function [knnModel_1, knnModel_2, knnModel_3] = createModel(path, txt_path, N, seed, K)
% Function to create the k-NN classification model for character detection from a text image. 
%   - reads the image and ASCII files, detects the characters and then finds 
%   the contours of each one. 
%   - then, creates a dataset for each type of character based on the numbers of
%   contours is has (1, 2 or 3), which contains a metric for each contour 
%   and the label ( ASCII code ) of the character.
%   - also, performs train-test split and trains a k-NN model for each dataset.
%   - then, it calculates the confusion matrix and weighted accuracy of each
%   model.
%
% Inputs:
%   path - a string containing the path to the text containing image
%   txt_path - a string containing the path to the .txt file containing the
%   text
%   N - a 1x3 array containing N1, N2, N3 values for number of points for
%   interpolation of contour describing points.
%   seed - an integer defining the seed for train-test split
%   K - a 1x3 array containing k1, k2, k3 values for each k-NN model.
%
% Outputs: 
%   knnModel_1 - the k-NN trained model for 1-contour character
%   classification
%   knnModel_2 - the k-NN trained model for 2-contour chars
%   knnModel_3 - the k-NN trained model for 3-contour chars



    img = imread(path);
    
    [m,n, numColor] = size(img);
    
    if numColor == 3
        gray_img = rgb2gray(img);
    else
        gray_img = img;
    end

    % Preprocessing 
    % If image is text1_v3 or text2_unrot, add a 10px white-padding in each side
    if strcmp(path, "text2_150dpi_unrot.png") || strcmp(path, "text1_v3.png")
        temp = uint8(zeros(m+20, n+20) + 255);
        temp(11:m+10, 11:n+10) = gray_img;
        gray_img = temp;
    end

    % Call splitWords function
    words = splitWords(gray_img);
    
    % Call splitCharacters function
    chars = splitCharacters(words);
    
    % Call createDataset function with N1, N2, N3 (number of points
    % describing each contour)
    [dataset_1, dataset_2, dataset_3] = createDataset(chars, txt_path, N);
    
    % Call splitDataset function
    % Define the random seed for train-test split for results repeatability
    [training_1, testing_1, training_2, testing_2, training_3, testing_3] = splitDataset(dataset_1, dataset_2, dataset_3, seed);
    
    % MODEL 1 ( characters with one contour )
    % Create X_train, y_train, X_test, y_test in appropriate forms for the
    X_train = training_1(:,1);
    y_train = uint8(cell2mat(training_1(:,2)));
    X_test = testing_1(:,1);
    y_test = uint8(cell2mat(testing_1(:,2)));
    
    train_m = size(X_train, 1);
    test_m = size(X_test, 1);
    
    train_n = length(cell2mat(X_train{1,1}));
    test_n = length(cell2mat(X_test{1,1}));
    
    X_train_app = zeros(train_m, train_n);
    X_test_app = zeros(test_m, test_n);
    
    for i = 1:train_m
        X_train_app(i, :) = cell2mat(X_train{i});
    end
    for i = 1:test_m
        X_test_app(i, :) = cell2mat(X_test{i});
    end
    
    % Create k-NN model
    % Define custom distance function using squared error
    squared_error_distance = @(X, Y) sum((X - Y).^2, 2);
    
    k = K(1);
    knnModel_1 = fitcknn(X_train_app, y_train, 'Distance', squared_error_distance, 'NumNeighbors', k);
    %knnModel = fitcknn(X_train_app, y_train);
    
    % Predict with the k-NN model
    y_pred = predict(knnModel_1, X_test_app);
    
    % Create char labels for visualization
    y_test_c = char(y_test);
    y_pred_c = char(y_pred);
    
    % Create and plot the confusion matrix of the model
    figure, clf;
    cm = confusionchart(y_test_c, y_pred_c);
    cm.Title = 'Model 1: Confusion Matrix';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    
    % Calculate the weighted accuracy of the model
    C = confusionmat(y_test_c, y_pred_c);
    
    class_counts = sum(C, 2);
    class_accuracies = diag(C) ./ class_counts;
    class_proportions = class_counts / sum(class_counts);
    weighted_accuracy1 = sum(class_accuracies .* class_proportions);
    fprintf('Model-1: one contour characters, N = %d, k = %d \n', N(1), k);
    fprintf('Weighted accuracy: %.2f\n', weighted_accuracy1);
    
    % MODEL 2 ( characters with two contours )
    % Create X_train, y_train, X_test, y_test in appropriate forms for the
    X_train = training_2(:,1);
    y_train = uint8(cell2mat(training_2(:,2)));
    X_test = testing_2(:,1);
    y_test = uint8(cell2mat(testing_2(:,2)));
    
    train_m = length(X_train);
    test_m = length(X_test);
    
    train_n = length(cell2mat(X_train{1,1}));
    test_n = length(cell2mat(X_test{1,1}));
    
    X_train_app = zeros(train_m, train_n);
    X_test_app = zeros(test_m, test_n);
    
    for i = 1:train_m
        X_train_app(i, :) = cell2mat(X_train{i, 1});
    end
    for i = 1:test_m
        X_test_app(i, :) = cell2mat(X_test{i, 1});
    end
    
    % Create k-NN model
    % Define custom distance function using squared error
    squared_error_distance = @(X, Y) sum((X - Y).^2, 2);
    
    k = K(2);
    knnModel_2 = fitcknn(X_train_app, y_train, 'Distance', squared_error_distance, 'NumNeighbors', k);
    
    % Predict with the k-NN model
    y_pred = predict(knnModel_2, X_test_app);
    
    % Create char labels for visualization
    y_test_c = char(y_test);
    y_pred_c = char(y_pred);
    
    % Create and plot the confusion matrix of the model
    figure, clf;
    cm = confusionchart(y_test_c, y_pred_c);
    cm.Title = 'Model 2: Confusion Matrix';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    
    % Calculate the weighted accuracy of the model
    C = confusionmat(y_test_c, y_pred_c);
    
    class_counts = sum(C, 2);
    class_accuracies = zeros(size(class_counts));  % Initialize accuracies
    for i = 1:numel(class_counts)
        if class_counts(i) ~= 0
            class_accuracies(i) = C(i,i) / class_counts(i);
        end
    end
    class_proportions = class_counts / sum(class_counts);
    weighted_accuracy2 = sum(class_accuracies .* class_proportions);
    fprintf('Model-2: two contour characters, N = %d, k = %d \n', N(2), k);
    fprintf('Weighted accuracy: %.2f\n', weighted_accuracy2);
    
    % MODEL 3 ( characters with three contours )
    % Create X_train, y_train, X_test, y_test in appropriate forms for the
    X_train = training_3(:,1);
    y_train = uint8(cell2mat(training_3(:,2)));
    X_test = testing_3(:,1);
    y_test = uint8(cell2mat(testing_3(:,2)));
    
    train_m = length(X_train);
    test_m = length(X_test);
    
    train_n = length(cell2mat(X_train{1,1}));
    test_n = length(cell2mat(X_test{1,1}));
    
    X_train_app = zeros(train_m, train_n);
    X_test_app = zeros(test_m, test_n);
    
    for i = 1:train_m
        X_train_app(i, :) = cell2mat(X_train{i, 1});
    end
    for i = 1:test_m
        X_test_app(i, :) = cell2mat(X_test{i, 1});
    end
    
    % Create k-NN model
    % Define custom distance function using squared error
    squared_error_distance = @(X, Y) sum((X - Y).^2, 2);
    
    k = K(3);
    knnModel_3 = fitcknn(X_train_app, y_train, 'Distance', squared_error_distance, 'NumNeighbors', k);
    
    % Predict with the k-NN model
    y_pred = predict(knnModel_3, X_test_app);
    
    % Create char labels for visualization
    y_test_c = char(y_test);
    y_pred_c = char(y_pred);
    
    % Create and plot the confusion matrix of the model
    figure, clf;
    cm = confusionchart(y_test_c, y_pred_c);
    cm.Title = 'Model 3: Confusion Matrix';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    
    % Calculate the weighted accuracy of the model
    C = confusionmat(y_test_c, y_pred_c);
    
    class_counts = sum(C, 2);
    class_accuracies = zeros(size(class_counts));  % Initialize accuracies
    for i = 1:numel(class_counts)
        if class_counts(i) ~= 0
            class_accuracies(i) = C(i,i) / class_counts(i);
        end
    end
    class_proportions = class_counts / sum(class_counts);
    weighted_accuracy3 = sum(class_accuracies .* class_proportions);
    fprintf('Model-3: three contours characters, N = %d, k = %d \n', N(3), k);
    fprintf('Weighted accuracy: %.2f\n', weighted_accuracy3);
    
    % Calculate total accuracy of model
    s1 = size(testing_1, 1);
    s2 = size(testing_2, 1);
    s3 = size(testing_3, 1);
    s = s1+s2+s3;
    
    total_acc = (s1/s)*weighted_accuracy1 + (s2/s)*weighted_accuracy2 + (s3/s)*weighted_accuracy3;
    fprintf('Total model%cs weighted accuracy: %.2f\n\n', "'",total_acc);
end