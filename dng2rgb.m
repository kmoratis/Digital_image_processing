function [Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, bayertype, method, M, N)
% Function implementing the .DNG to RGB transform, by calling the functions 
% wbmask() for white-balance handling, demosaicing_nearest() or demosaicing_linear() according to the method arg
% and space_transform() for space transformations.
%
% Inputs:
%   rawim - the input raw image in Bayer pattern
%   XYZ2Cam - 3x3 matrice for XYZ_to_Cam transform
%   wbcoeffs - 1x3 vector for white-balance handling
%   bayertype - a string specifying the Bayer pattern ('RGGB', 'BGGR', 'GBRG', 'GRBG')
%   method - a string specifying the interpolation method ('nearest', 'linear')
%   M, N - the size of the output R, G, B matrices
%
% Outputs:
%   Csrgb - The output matrice in R, G, B form
%   Clinear - The output matrice in Rlinear, Glinear, Blinear form
%   Cxyz - The output matrice in X, Y, Z form
%   Ccam - The output matrice in Rcam, Gcam, Bcam form

    % White balancing using wbmask
    mask = wbmask(size(rawim,1), size(rawim,2), wbcoeffs, bayertype);
    balanced_bayer = rawim .* mask;

    % Demosaicing the Bayer-type raw image to RGB 
    if strcmp(method, 'nearest')
        Ccam = demosaicing_nearest(balanced_bayer, bayertype, M, N);
    elseif strcmp(method, 'linear')
        Ccam = demosaicing_linear(balanced_bayer, bayertype, M, N);
    else
        error('Invalid interpolation method type');
    end
    figure('Name','Ccam');
    imshow(Ccam);    
    title('Image 1');

    % space transform
    [Cxyz, Clinear, Csrgb] = space_transform(Ccam, XYZ2Cam);
    figure('Name','Cxyz');
    imshow(Cxyz);
    title('Image 2');
    figure('Name','Clinear');
    imshow(Clinear);
    title('Image 3');
    figure('Name','Csrgb');
    imshow(Csrgb);
    title('Image 4');
end

