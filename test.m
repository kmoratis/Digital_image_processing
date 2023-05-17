% Create dataset test
close all;

%path = "text2_150dpi_unrot.png";
path = "text1_v3.png";
opening_rad = 3;

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

% If image is text2_150dpi_rot, make edges (2px each) white
if strcmp(path, "text2_150dpi_rot.png")
    
    % Preprocessing: If image is text2, make edges (2px each) white
    gray_img(1:2, :) = 255;
    gray_img(m-1:m, :) = 255;
    gray_img(:, 1:2) = 255;
    gray_img(:, n-1:n) = 255;
end

% Call splitWords function
words = splitWords(gray_img);

% Call splitCharacters function
chars = splitCharacters(words);

% Call createDataset function
[dataset_1, dataset_2, dataset_3] = createDataset(chars);

% Call splitDataset function
[training_1, testing_1, training_2, testing_2, training_3, testing_3] = splitDataset(dataset_1, dataset_2, dataset_3);

% Create X_train, y_train, X_test, y_test in appropriate forms for the
% model 1
X_train = training_1(:,1);
y_train = uint8(cell2mat(training_1(:,2)));
X_test = testing_1(:,1);
y_test = uint8(cell2mat(testing_1(:,2)));

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
k = 3;
%knnModel = fitcknn(X_train, y_train, 'Distance', 'mse', 'NumNeighbors', k, 'Standardize', 1);
knnModel = fitcknn(X_train_app, y_train);

% Train the k-NN model
knnModel = train(knnModel, X_train_app, y_train);

% Predict with the k-NN model
y_pred = predict(knnModel, X_test_app);


