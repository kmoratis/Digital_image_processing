function colormask = wbmask(m, n, wbmults, align)
% Implementing a white-balance multiplicative mask for an m-by-n image 
%
% Inputs: 
%   m, n - the image size
%   wbmults - RGB white balance multipliers wbmults = [R_scale, G_scale, B_scale]
%   align - a sting indicating Bayer arrangement: 'RGGB', 'BGGR', 'GRBG', 'GBRG'
%
% Outputs:
%   colormask - the desired mask

    colormask = wbmults(2)*ones(m,n); %Initialize to all green values
    switch align
        case 'RGGB'
            colormask(1:2:end,1:2:end) = wbmults(1);    %r
            colormask(2:2:end,2:2:end) = wbmults(3);    %b
        case 'BGGR'
            colormask(2:2:end,2:2:end) = wbmults(1);    %r
            colormask(1:2:end,1:2:end) = wbmults(3);    %b
        case 'GRBG'
            colormask(1:2:end,2:2:end) = wbmults(1);    %r
            colormask(1:2:end,2:2:end) = wbmults(3);    %b
        case 'GBRG'
            colormask(2:2:end,1:2:end) = wbmults(1);    %r
            colormask(1:2:end,2:2:end) = wbmults(3);    %b
        otherwise
            disp('Invalid Bayer pattern type');
    end
end