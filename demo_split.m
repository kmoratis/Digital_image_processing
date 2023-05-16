%clear;
%clc;
%close all;

gray_img = imread('text1_unrot.png');

% Call splitCharacters function
words = splitWords(gray_img);

%%%%%% test below

total_chars_found = 0;
all_chars = {};

for i = 1:length(words)-1
    curr_word = cell2mat(words(i)); 

    % Pad with white 1 px each side
    [m, n] = size(curr_word);
    temp = uint8(zeros(m+2, n+2) + 255);
    temp(2:m+1, 2:n+1) = curr_word(:, :);
    curr_word = temp;

    % Find projection to the horizontal axis of the word's pixels
    horizontal_proj = sum(curr_word, 1);

    % Get white column projection and create a threshold value
    white_proj = horizontal_proj(1);
    threshold = 0.953 * white_proj; %0.946 working pretty good

    start = 0;
    stop = 0;
    % Detect word's chars
    for j = 1:length(horizontal_proj)
        % Get current chars projection
        char_proj = horizontal_proj(1, j);

        % Found start
        if char_proj <= threshold && start==0
            start = j;
            total_chars_found = total_chars_found + 1;
        % Found stop 
        elseif char_proj >= threshold && start~=0
            stop = j;
            char_img = curr_word(:, start:stop);
            all_chars(1, total_chars_found) = {char_img};
            start = 0;
            stop = 0;
        end

    end
end

offset = 500;
figure(), clf;
for i = 1:20
    subplot(1, 20, i);
    imshow(cell2mat(all_chars(1, offset+i)));
end