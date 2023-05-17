function [dataset_1, dataset_2, dataset_3] = createDataset(chars, N)
%CREATEDATASET Summary of this function goes here
%   Detailed explanation goes here

    number_chars = length(chars);

    txt_chars = fileread('text1_v3.txt', 'Encoding','UTF-8');
    num_txt_chars = length(txt_chars);
    labels = zeros(1, num_txt_chars);
    labels_found = 0;

    % Extract only letters and numbers, '.', ',', '-', '!', '?', ';'
    cleaned_txt = regexp(txt_chars, '[a-zA-Z0-9.,-!?;]', 'match');
    cleaned_txt = [cleaned_txt{:}];

    num_cleaned = length(cleaned_txt);
    
    % Check if character detect was correct
    if num_cleaned ~= number_chars
        disp('Something went wrong with character split: ');
        disp(['Characters found: ', num2str(number_chars)]);
        disp(['Actual characters: ', num2str(num_cleaned)]);
    end

    % Calculate the ASCII values
    ascii_values = double(cleaned_txt);

    % Create three datasets to store each character type (1-2-3 contours)
    dataset_1 = cell(number_chars, 2);
    found_1 = 0;
    dataset_2 = cell(number_chars, 2);
    found_2 = 0;
    dataset_3 = cell(number_chars, 2);
    found_3 = 0;

    for i = 1:number_chars
        c = cell2mat(chars(i));

        % call getcontour function
        char_contours = getcontour(c, false);

        % call calculateDescriptor function, with N1=N2=N3 = 500
        N = [500 500 500];
        char_descriptors = calculateDescriptor(char_contours, N);
        
        % character with one contour
        if length(char_contours) == 1
            found_1 = found_1 + 1;
            dataset_1(found_1, 1) = {char_descriptors};
            dataset_1(found_1, 2) = {ascii_values(i)};
        
        % character with two contours
        elseif length(char_contours) == 2
            found_2 = found_2 + 1;
            dataset_2(found_2, 1) = {char_descriptors};
            dataset_2(found_2, 2) = {ascii_values(i)};
        
        % character with three contours
        elseif length(char_contours) == 3
            found_3 = found_3 + 1;
            dataset_3(found_3, 1) = {char_descriptors};
            dataset_3(found_3, 2) = {ascii_values(i)};
        else
            disp('Something went wrong');
        end     
    end

    % keep the usefull part of the dataset arrays
    dataset_1 = dataset_1(1:found_1, :);
    dataset_2 = dataset_2(1:found_2, :);
    dataset_3 = dataset_3(1:found_3, :);

end