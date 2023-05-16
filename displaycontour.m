function out = displaycontour(c, m, n)
% Function to display the contours of the letter found, each one is
% separate color.
%
% Inputs:
%   c - The 3 cell array containing one Nx2 array for each contour found,
%   else it contains 0
%   m,n - The size of the initial image
%
% Output:
%   out - the image containing the contours, each in different color

    % Create RGB image sized m+2, n+2 (original with 1 px pad)
    m = m+2;
    n = n+2;
    out = uint8(zeros(m, n, 3) + 255);

    [~, w] = size(c); % w should b 3, just for more generic sol

    for i = 1:w
        if i==1
            color = [0 0 255]; % Outside contour blue
        elseif i==2
            color = [255 0 0]; % Second inward red
        elseif i==3
            color = [0 255 0]; % Third (inside) green
        else 
            disp('Something went wrong');
        end

        % If all contours displayed, the cell contains only 0
        ci = cell2mat(c(i));
        if ci == 0
            break;
        end

        % Get number of points describing the contour
        [N, ~] = size(ci);
        for j = 1:N
            coo = ci(j, :);
            out(coo(1), coo(2), :) = color(1,:);
        end
    end

    % Create figure and display the image with different contours colors
    figure(), clf;
    imshow(out, 'InitialMagnification', 'fit');
    title('Character contours')

end