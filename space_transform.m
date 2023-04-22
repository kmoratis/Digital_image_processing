function [Cxyz, Clinear, Csrgb] = space_transform(Ccam, XYZ2Cam)
% Performs space transformations from Ccam to Cxyz, Clinear, Csrgb formats.
%  
% Inputs:
%   Ccam - matrice containing R, G, B channels.
%   XYZ2Cam - 3x3 matrice for XYZ_to_Cam transform
%
% Outputs:
%   Cxyz - matrice containing X, Y, Z channels
%   Clinear - matrice containing Rlinear, Glinear, Blinear channels
%   Csrgb - matrice containing R, G, B channels
    
    % Cam to XYZ transform
    Cam2XYZ = XYZ2Cam^-1;
    Cxyz = apply_cmatrix(Ccam, Cam2XYZ);
    Cxyz = max(0,min(Cxyz,1)); % Clip image b/w 0-1

    % Cam to linear transform
    % Define XYZ2RGB array
    XYZ2RGB = [3.2406 -1.5372 -0.4986; -0.9689 1.8758 0.0415; 0.0557 -0.2040 1.0570];
    RGB2XYZ = XYZ2RGB^-1;
    RGB2Cam = XYZ2Cam * RGB2XYZ;
    RGB2Cam = RGB2Cam ./ repmat(sum(RGB2Cam,2),1,3); % Normalize rows to 1
    Cam2RGB = RGB2Cam^-1;
    Clinear = apply_cmatrix(Ccam, Cam2RGB);
    Clinear = max(0, min(Clinear, 1)); % Clip image b/w 0-1

    % linear to sRGB transform
    % brightness and gamma correction
    grayim = rgb2gray(Clinear);
    grayscale = 0.25/mean(grayim(:));
    bright_srgb = min(1, Clinear*grayscale);

    Csrgb = bright_srgb.^(1/2.2);

end

