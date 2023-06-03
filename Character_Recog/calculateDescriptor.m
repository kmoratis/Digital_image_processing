function char_descriptors = calculateDescriptor(char_contours, N)
% Function calculating the descriptors for each contour of a character
%
% Inputs:
%   char_contours - A cell array containing the contour points of each
%   contour of the character
%   N - a 1x3 array, containing N1, N2, N3 values ( number of points
%   describing each contour ), for interpolation
%
% Outputs:
%   char_descriptors - A cell array containing one cell for each contour
%   found. Each cell contains N1, N2 or N3 values, which is the descriptor
%   of the character's contour.


    num_contours = length(char_contours);
    char_descriptors = cell(1, num_contours);

    for i = 1:num_contours
        current_char = cell2mat(char_contours(i));

        % find the complex sequence
        x = current_char(:, 1);
        y = current_char(:, 2);
        n = length(x);

        % interp1 needs at least 2 points in each dimension
        if n > 1
            % desired number of points for the descriptor (different for each
            % contour)
            number_of_points = N(i);
    
            query_points = linspace(1, n, number_of_points);
    
            interp_x = interp1(1:n, x, query_points);
    
            interp_y = interp1(1:n, y, query_points);
            
            interp_r = interp_x + 1i*interp_y; %1i is j ( imaginary unit )

            % calculate DFT using FFT
            R = abs(fft(interp_r));

            % exclude the first value of DFT
            R = R(2:length(R));
    
            char_descriptors(i) = {R};

        % if n == 1
        else 
            interp_r = x + 1i*y; % not actually interpolated, but keeping the symbolism 

            char_descriptors(i) = {interp_r};

            % calculate DFT using FFT, and then keep the absolute
            R = abs(fft(interp_r));

            % exclude the first value of DFT
            R = R(2:length(R));
    
            char_descriptors(i) = {R};
        end
    end
end
