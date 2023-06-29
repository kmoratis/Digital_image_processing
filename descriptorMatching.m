function matchingPoints = descriptorMatching(im1, im2, points1, points2, percentageThreshold, desc_upgrade, imshow)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    figure, clf;
    imshow(im2);

    % Descriptor parameters
    rhom = 2; % Min rad
    rhoM = 12; % Max rad
    rhostep = 1; % Rad step
    N = 40; % Number of points
    
    % Calc descriptor vector size
    s = 1 + (rhoM - rhom) / rhostep;

    % Create array to store actual corner points ( descriptor not empty )
    cor_1 = zeros(size(points1));
    cor1_found = 0;
    
    % Create array to store im1 descriptors
    N1 = size(points1, 1);
    desc1 = zeros(s, N1);
    
    % Calculate im1 corners descriptors
    for i = 1:N1
        p = points1(i, :);
    
        if desc_upgrade == true % use desc_upgrade verion of the descriptor
            d = myLocalDescriptorUpgrade(im1, p, rhom, rhoM, rhostep, N);
            if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
                cor1_found = cor1_found + 1;
                desc1(:, cor1_found) = d(:);
                cor_1(cor1_found, :) = p;
            end
        else
            d = myLocalDescriptor(im1, p, rhom, rhoM, rhostep, N);
            if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
                cor1_found = cor1_found + 1;
                desc1(:, cor1_found) = d(:);
                cor_1(cor1_found, :) = p;
            end
        end
    end
    
    % keep the usefull part of the array
    cor_1 = cor_1(1:cor1_found, :);
    desc1 = desc1(:, 1:cor1_found);
    
    % Create array to store actual corner points ( descriptor not empty )
    cor_2 = zeros(size(points2));
    cor2_found = 0;
    
    % Create array to store im2 descriptors
    N2 = size(points2, 1);
    desc2 = zeros(s, N2);
    
    % Calculate im2 corners descriptors
    for i = 1:N2
        p = points2(i, :);
    
        if desc_upgrade == true % use desc_upgrade verion of the descriptor
            d = myLocalDescriptorUpgrade(im2, p, rhom, rhoM, rhostep, N);
            if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
                cor2_found = cor2_found + 1;
                desc2(:, cor2_found) = d(:);
                cor_2(cor2_found, :) = p(1, :);
            end
        else
            d = myLocalDescriptor(im2, p, rhom, rhoM, rhostep, N);
            if ~isempty(d) && sum(d)~=0 % if descriptor is not empty vector
                cor2_found = cor2_found + 1;
                desc2(:, cor2_found) = d(:);
                cor_2(cor2_found, :) = p;
            end
        end
    end
    
    % keep the usefull part of the array
    cor_2 = cor_2(1:cor2_found, :);
    desc2 = desc2(:, 1:cor2_found);
    
    %% Descriptor Matching
    N1 = size(desc1, 2);
    N2 = size(desc2, 2);
    
    desc_match = zeros(N1, N2);
    
    % Create N1xN2 array containing the descriptor simillarity of each point
    % pair
    for i = 1:N1
        nm = norm(desc1(:,i));
        for j = 1:N2
            desc_match(i, j) = norm(desc1(:, i) - desc2(:, j)) / nm; % normalised distance
        end
    end
    
    % Create matchingPoints array, containing the n pairs that have a
    % simillarity above percentageThreshold value
    
    % Flatten the descriptor matching array and get the indices
    [flat_array, indices] = sort(desc_match(:), 'ascend'); % Looking for the min distance
    
    % Find the number of elements, representing the percentageThreshold % of the total values
    num_elements = ceil(length(flat_array) * percentageThreshold);
    
    % Select the largest elements and their indices
    sm_el = flat_array(1:num_elements);
    sm_in = indices(1:num_elements);
    
    % Convert the linear indices back to row and column ones
    [row_i, col_i] = ind2sub(size(desc_match), sm_in);
    
    matchingPoints = [row_i col_i];
    
    % Find best match and plot it 
    p1_idx = matchingPoints(20,1); % best_i -> index of best point of img1
    p2_idx = matchingPoints(20,2); % best_j
    
    p1 = cor_1(p1_idx, :);
    p2 = cor_2(p2_idx, :);
    
    if imshow
        % Display the corners on the original image
        figure, clf;
        subplot(1,2,1);
        imshow(im1);
        hold on
        plot(p1(1, 2), p1(1, 1), 'r+', 'LineWidth', 10);
        hold off
        
        subplot(1,2,2);
        imshow(im2);
        hold on
        plot(p2(1,2), p2(1,1), 'r+', 'LineWidth', 10);
        hold off
    end

end