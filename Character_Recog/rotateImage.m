function y = rotateImage(x, angle)
% Rotate the input grayscale image x by the specified angle (in degrees) counter-clockwise, with white
% padding 
%
% Inputs:
%   x - the grayscale image to be rotated
%   angle - the angle of rotation in degrees
%   interp_method - the interpolation method for the imrotate
%
% Output:
%   y - the rotated grayscale image with white padding

    switch mod(angle, 360)
        % Special cases
        case 0
            y = x;
        case 90
            y = rot90(x);
        case 180
            y = x(end:-1:1, end:-1:1);
        case 270
            y = rot90(x(end:-1:1, end:-1:1));
    
        % General rotations
        otherwise

            % Convert to radians and create transformation matrix R
            rads = angle*pi/180;
            R = [+cos(rads) +sin(rads); -sin(rads) +cos(rads)];
    
            % Figure out the size of the original image
            [m,n,p] = size(x);
            % Calculate corners of the rotated image
            dest = round( [1 1; 1 n; m 1; m n]*R );
            % Adjust the coorfinates of the corners to ensure that all
            % values are at least 1
            dest = bsxfun(@minus, dest, min(dest)) + 1;
            % Create a white image for rotated image
            y = zeros([max(dest) p],class(x)) + 255;
    
            % Map all pixels of the transformed image to the original image
            for ii = 1:size(y,1)
                for jj = 1:size(y,2)
                    source = ([ii jj]-dest(1,:))*R.';
                    if all(source >= 1) && all(source <= [m n])
    
                        % Get all 4 surrounding pixels
                        C = ceil(source);
                        F = floor(source);
    
                        % Compute the relative areas
                        A = [...
                            ((C(2)-source(2))*(C(1)-source(1))),...
                            ((source(2)-F(2))*(source(1)-F(1)));
                            ((C(2)-source(2))*(source(1)-F(1))),...
                            ((source(2)-F(2))*(C(1)-source(1)))];
    
                        % Extract value and re-scale it relative to area
                        cols = bsxfun(@times, A, double(x(F(1):C(1),F(2):C(2),:)));
    
                        % Assign                     
                        y(ii,jj,:) = sum(sum(cols),2);
    
                    end
                end
            end        
    end

end
