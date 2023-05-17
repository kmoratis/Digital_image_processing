function word_imgs = splitWords(img)
% This function takes a grayscale text image as an input ( after rotation reverse, if needed )
% and detects all its words. Firstly, it detects all lines, and then for
% each line it detects its words.
% Note: If two or more words are underlined as a sentence, meaning
% there is not a row of pixels between them that is all white, then the algorithm
% detects them as one word.
%
% Input:
% img - the grayscale text image, after rotation reverse
%
% Outputs:
% word_imgs - a cell array, containing each word's image


    % Line detection using the projection to the vertical axis
    % Each line will have a white padding on top and bottom
 
    % Calculate the vertical projection
    vertical_proj = sum(img, 2);
    
    % Get white line vertical projection ( we ensured that the first row will be
    % all whites in the previous step )
    white_line = vertical_proj(1);
    
    % Create zero array to store line starting and ending points ( we suppose
    % that the lines will be less that 100 )
    lines = zeros(100, 2);
    lines_found = 0;
    start = 0;
    ending = 0;

    % Create cell array to store found words
    word_imgs = {};
    
    % For each line (not-white), save its starting and ending point
    for i = 2:length(vertical_proj)
        line_proj = vertical_proj(i);
       
        % Found start
        if (line_proj ~= white_line) && start==0 
            lines_found = lines_found+1;
            start = i;
        end
        % Found ending
        if (line_proj == white_line) && start~=0
            ending = i-1;
            lines(lines_found, :) = [start, ending];
            start = 0;
            ending = 0;
        end
    end
    
    % Keep the usefull part of the array
    lines = lines(1:lines_found, :);
    %disp(['Lines found: ', num2str(lines_found)]);

    total_words_found = 0;
    
    % Finding an approximation of the space length in columns, so the
    % funtion is as responsive as possible (still not perfect)
    [~, n] = size(img);
    space_length = ceil(n/110);

    % Detect the words of each line, using the projection to the horizontal axis
    for j = 1:lines_found
        start = lines(j, 1);
        stop = lines(j, 2);
     
        % Calculate horizontal projection
        curr_line = img(start:stop, :);
        horizontal_proj = sum(curr_line, 1);

        % Get white pixels row projection
        white_col = horizontal_proj(1);
     
        % Create zero array to store word starting and ending points ( we suppose
        % that the words in a line will be less that 100 )
        w_start = 0;
        w_ending = 0;
        
        % For each word save an image containing it
        for i = 2:length(horizontal_proj)
        
            col_proj = horizontal_proj(i);
            % Found start
            if (col_proj ~= white_col) && w_start==0 
                total_words_found = total_words_found + 1;
                w_start = i-1;
            end
            % Found ending - at least (space_length) white rows
            if (col_proj == white_col) && w_start~=0
                found_space = 1;

                % Check next rows
                % Space_length || end of line
                resp_length = min(space_length, n-i);
                for k = 1:resp_length
                    if (horizontal_proj(i+k) == white_col)
                        continue;
                    else
                        found_space = 0;
                        break;
                    end
                end
                % If all columns after first white are also white (found_space),
                % then add word to found words
                if found_space==1
                    w_ending = i;
                    w_img = curr_line(:, w_start:w_ending);
                    word_imgs(total_words_found) = {w_img};

                    w_start = 0;
                    w_ending = 0;
                    found_space = 0;
                    i = i + resp_length; % skip next iterations
                end
            end
        end     
    end
end