function readText(img_path, txt_path, knnModel_1, knnModel_2, knnModel_3, N)
% Function for detecting characters and predicting their class with the
% given the trained knnModels.
% Note: image rotation takes quite some time, so you can comment it for
% faster testing ( and provide unrotated image instead )
% 
% Inputs:
%   img_path - string containing the path to the text image
%   txt_path - string containing the path to the text file
%   knnModel_1 - trained k-NN model for 1-contour character classification
%   knnModel_2 - trained k-NN model for 2-contour character classification
%   knnModel_3 - trained k-NN model for 3-contour character classification
%

    % Read image and transform it to grayscale if needed
    img = imread(img_path);
    
    [~,~, numColor] = size(img);
    
    if numColor == 3
        gray_img = rgb2gray(img);
    else
        gray_img = img;
    end

    % Preprocessing
    % Only for text1_rot to remove black padding
    if strcmp(img_path, "text2_150dpi_rot.png")
        % Preprocessing to remove the black padding
        [m, n] = size(gray_img);
        gray_img(1:9, :) = 255; % top
        gray_img(m-9:m, :) = 255; % bottom
        gray_img(:, 1:7) = 255; % left 
        gray_img(:, n-7:n) = 255; % right
    end
    
    % Only for text1_v3
    if strcmp(img_path, "text1_v3.png")
        % Add one pixel white-padding to each side
        temp = uint8(zeros(rows+2, cols+2) + 255);
        temp(2:rows+1, 2:cols+1) = gray_img(:,:);
        gray_img = temp;
        %imshow(gray_img);
    end

    % Find the rotation angle of the input image
    precision = 15; % number of points to be examined in linear search
    fig_show = false;
    angle = findRotationAngle(gray_img, precision, fig_show);

    % Call rotateImage function to reverse the image rotation
    unrot_img = rotateImage(gray_img, -angle);
    
    % Call splitWords function
    words = splitWords(unrot_img);
    
    % Call splitCharacters function
    chars = splitCharacters(words);
    
    % Call createDataset function with N1, N2, N3 (number of points
    % describing each contour)
    [testset_1, testset_2, testset_3] = createDataset(chars, txt_path, N);

    % Evaluate model 1
    X_test = testset_1(:,1);
    y_test = uint8(cell2mat(testset_1(:,2)));
    
    test_m = size(X_test, 1);
    test_n = length(cell2mat(X_test{1,1}));
    
    X_test_app = zeros(test_m, test_n);
    for i = 1:test_m
        X_test_app(i, :) = cell2mat(X_test{i});
    end
    
    % Predict with the k-NN model
    y_pred = predict(knnModel_1, X_test_app);
    
    % Create char labels for visualization
    y_test_c = char(y_test);
    y_pred_c = char(y_pred);
    
    % Create and plot the confusion matrix of the model
    figure, clf;
    cm = confusionchart(y_test_c, y_pred_c);
    cm.Title = 'Model 1: Confusion Matrix (Test set (text 2)';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    
    % Calculate the weighted accuracy of the model
    C = confusionmat(y_test_c, y_pred_c);
    
    fprintf('Evaluating Models using test set: (%s)\n', img_path);
    class_counts = sum(C, 2);
    class_accuracies = zeros(size(class_counts));  % Initialize accuracies
    for i = 1:numel(class_counts)
        if class_counts(i) ~= 0
            class_accuracies(i) = C(i,i) / class_counts(i);
        end
    end
    class_proportions = class_counts / sum(class_counts);
    weighted_accuracy1 = sum(class_accuracies .* class_proportions);
    fprintf('Model-1: Weighted accuracy: %.2f\n', weighted_accuracy1);
    
    % Evaluate model 2
    X_test = testset_2(:,1);
    y_test = uint8(cell2mat(testset_2(:,2)));
    
    test_m = size(X_test, 1);
    test_n = length(cell2mat(X_test{1,1}));
    
    X_test_app = zeros(test_m, test_n);
    for i = 1:test_m
        X_test_app(i, :) = cell2mat(X_test{i});
    end
    
    % Predict with the k-NN model
    y_pred = predict(knnModel_2, X_test_app);
    
    % Create char labels for visualization
    y_test_c = char(y_test);
    y_pred_c = char(y_pred);
    
    % Create and plot the confusion matrix of the model
    figure, clf;
    cm = confusionchart(y_test_c, y_pred_c);
    cm.Title = 'Model 2: Confusion Matrix (Test set (text 2)';
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
    fprintf('Model-2: Weighted accuracy: %.2f\n', weighted_accuracy2);
    
    % Evaluate model 3
    X_test = testset_3(:,1);
    y_test = uint8(cell2mat(testset_3(:,2)));
    
    test_m = size(X_test, 1);
    test_n = length(cell2mat(X_test{1,1}));
    
    X_test_app = zeros(test_m, test_n);
    for i = 1:test_m
        X_test_app(i, :) = cell2mat(X_test{i});
    end
    
    % Predict with the k-NN model
    y_pred = predict(knnModel_3, X_test_app);
    
    % Create char labels for visualization
    y_test_c = char(y_test);
    y_pred_c = char(y_pred);
    
    % Create and plot the confusion matrix of the model
    figure, clf;
    cm = confusionchart(y_test_c, y_pred_c);
    cm.Title = 'Model 3: Confusion Matrix (Test set (text 2)';
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
    fprintf('Model-3: Weighted accuracy: %.2f\n', weighted_accuracy3);
    
    % Calculate total accuracy of model
    s1 = size(testset_1, 1);
    s2 = size(testset_2, 1);
    s3 = size(testset_3, 1);
    s = s1+s2+s3;
    
    total_acc = (s1/s)*weighted_accuracy1 + (s2/s)*weighted_accuracy2 + (s3/s)*weighted_accuracy3;
    fprintf('Total model%cs weighted accuracy: %.2f\n\n', "'",total_acc);

end