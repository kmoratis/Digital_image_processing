function matchingPoints = descriptorMatching(desc1, desc2, percentageThreshold)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
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

    matchingPoints = matchingPoints'; % 2xn array, containing n pairs of matched points ( ascending sorted )

end