function matchingPoints = descriptorMatching(desc1, desc2, percentageThreshold)
% Function for finding the salient point pairs with the best matches,
% according to the normalized Euclidean distance of theirs descriptor
% vectors.
%
% Inputs:
%   desc1 - image_1 salient points descriptors (N1xN2 array, where N1 is
%       the length of each pixel's descriptor, and N2 is the number of salient
%       points)
%   desc2 - image_2 salient points descriptors (N1xN2) array
%   percentageThreshold - a float number, for keeping only the
%       (percentageThreshold*100) % best matches, from the total ones.
%
% Output:
%   matchingPoints - a 2xn array of pixel pairs that are most likely to
%       be matched.
    
    % Descriptor Matching
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