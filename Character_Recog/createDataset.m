function [dataset_1, dataset_2, dataset_3] = createDataset(chars, txt_path, N)
% Function to create the three datasets containing the character images
%
% Inputs:
%   chars - A cell array. Each cell contains a grayscale image of a single
%   character
%   txt_path - a string containing the path to the .txt file with the input
%   text
%   N - a 1x3 array, containing N1, N2, N3 values for each contour ( how
%   many points describes each contour ), for interpolation
%
% Outputs:
%   dataset_1, dataset_2, dataset_3 - three cell arrays with sizes M1x2,
%   M2x2, M3x2, where M1, M2, M3 is the number of characters found of each
%   category ( 1, 2 or 3 contours ). The first column contains a cell
%   array, containing the descriptors of each contour. The second column
%   contains an ASCII code, indicating the label of the data point.

    number_chars = length(chars);

    txt_chars = fileread(txt_path, 'Encoding', 'UTF-8');
    num_txt_chars = length(txt_chars);
    labels = zeros(1, num_txt_chars);
    labels_found = 0;

    % Extract only letters and numbers, '.', ',', '-', '!', '?', ';', '(',
    % ')', ' ' '
    cleaned_txt = regexp(txt_chars, "[a-zA-Z0-9().,\-!?;()\']", 'match');
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

        % call calculateDescriptor function
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
