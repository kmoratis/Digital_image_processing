function c = getcontour(x, show_img)
% Given an input image (grayscale) containing a character, the function 
% performs the appropriate morphological transforms to the image,
% and then it finds all the contours of the character.
%
% Inputs:
%   x - A grayscale image containing a character
%   show_img - A boolean variable, to handle whether the figures will be shown
%
% Output:
%   c - A cell array containing an array of Nx2 points for each contour 
%   found, where N is the number of points describing the contour

    % Convert to binary
    binary_img = imbinarize(x);
    
    % Convert to dilated
    dilated_img = imdilate(binary_img, strel('disk', 1));
    
    % Remove the binary image from the dilated one
    contour_img = dilated_img - binary_img;
    
    % Apply thinning to make border size equal to 1px.
    thinned_img = bwmorph(contour_img, 'thin', Inf);

    % Display input img and transformed one
    if show_img
        figure(), clf;
        sgtitle('Image preprocessing for contour find') 
        subplot(1, 2, 1);
        imshow(x);
        title('Input image');
        subplot(1, 2, 2);
        imshow(thinned_img);
        title({'Transformed image', 'and Thinned to border size=1px'});
    end

    % Find thinned image size
    [m, n] = size(thinned_img);
    
    % Give 1px black-padding(same as background) to thinned image to ensure that 
    % the contour of the character in not at the edges
    temp = zeros(m+2, n+2);
    temp(2:m+1, 2:n+1) = thinned_img;
    thinned_img = temp;
         
    % Find new size
    [m, n] = size(thinned_img);
    
    % Make a copy of the thinned img to change its content
    thincp = thinned_img;
    
    % Display transformed image to new figure
    if show_img
        figure(), clf;
        subplot(2, 3, 2);
        sgtitle('Contour finding process') 
        imshow(thincp);
        title('Original Image');
    end
    
    % Initialize variables
    contour1 = 0;
    contour2 = 0;
    contour3 = 0;
    contours_found = 0;

    % Each iteration, finds one contour. Inward direction.
    % Do until all found.
    while(true)
    
        % When all contours found, exit
        if find(thincp) 
            contours_found = contours_found + 1;
        else
            break;
        end
    
        % Find boundary top-left pixel, as a starting point for the algorithm
        [bound_rows, bound_cols] = find(thincp);
        
        min_row = min(bound_rows);
        topleft_pixel_index = find(bound_rows == min_row, 1);
        
        min_col = bound_cols(topleft_pixel_index);
        
        % Create empty array for contour points with size N == m*n/2. ( Could possibly be picked even smaller ) 
        N = round(m*n/2);
        contour = zeros(N, 2);
        
        % Insert starting point of the contour, and delete it in the copy of the
        % thin image
        contour(1, :) = [min_row, min_col];
        thincp(min_row, min_col) = 0;
        
        % Find next point, with respect to right - right,bottom - bottom - left,bottom
        if thincp(min_row, min_col+1) == 1 % Check right neighbour
            contour(2, :) = [min_row, min_col+1];
            thincp(min_row, min_col+1) = 0;
            current = [min_row, min_col+1];
        elseif thincp(min_row+1, min_col+1) == 1 % Check bottom-right neighbour
            contour(2, :) = [min_row+1, min_col+1];
            thincp(min_row+1, min_col+1) = 0;
            current = [min_row+1, min_col+1];
        elseif thincp(min_row+1, min_col) == 1 % Check bottom neighbour
            contour(2, :) = [min_row+1, min_col];
            thincp(min_row+1, min_col) = 0;
            current = [min_row+1, min_col];
        elseif thincp(min_row+1, min_col-1) == 1 % CHeck bottom-left neighbour
            contour(2, :) = [min_row+1, min_col-1];
            thincp(min_row+1, min_col-1) = 0;
            current = [min_row+1, min_col-1];
        else
            disp('Something went wrong, on calculating first two elements');
        end
           
        point_num = 2;
        
        % Iterate through the next points of the contour, until no more neighbours
        % (Full Loop)
        while(true)
            % Find coordinates of neighbouring points 
            neighbors_idx = [current(1)-1 current(2)-1;
                            current(1)-1 current(2);
                            current(1)-1 current(2)+1;
                            current(1) current(2)-1;
                            current(1) current(2)+1;
                            current(1)+1 current(2)-1;
                            current(1)+1 current(2);
                            current(1)+1 current(2)+1];
    
            found_neigh = 0;
        
            % Find which one of them is 1 and add it to the contour array
            for i = 1:size(neighbors_idx, 1)
        
                coo = neighbors_idx(i, :);
                if thincp(coo(1), coo(2)) == 1
                    % Add point to contour array
                    contour(point_num+1, :) = neighbors_idx(i,:);
                    % Make point black
                    thincp(coo(1), coo(2)) = 0;
                    % Make it the current point
                    current = neighbors_idx(i, :);
                    % Increase number of points by 1
                    point_num = point_num + 1;
                    % Change found to 1
                    found_neigh = 1;
                    break;
                end
            end
        
            % If there are no more neighbouring points, the iteration was done,
            % exit
            if found_neigh == 0
                % Keep only the usefull part of the contour array
                contour = contour(1:point_num, :);
                break;
            end
        end
    
        % Display current version thinned image (without already found contours) 
        if show_img
            subplot(2, 3, 3 + contours_found);
            imshow(thincp);
            title(['After iteration ', num2str(contours_found)]);
        end
    
        % Store the contour to the appropriate position, inward direction.
        if contours_found == 1
            contour1 = contour;
        elseif contours_found == 2
            contour2 = contour;
        elseif contours_found == 3
            contour3 = contour;
        else
            disp('Something went wrong.');
        end
    end
    % Return the cell array containing each contour found
    c = {contour1, contour2, contour3};

end

