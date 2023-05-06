function [rawimg, XYZ2Cam, wbcoeffs] = readdng(filename)
% Reads RAW image (.DNG) and returns the useful part of the sensor value array, with white and black levels set.
%
% Inputs:
%   filename - a string containing the name (or path) of a .DNG file
%
% Outputs:
%   rawim - an array containing the sensor measurements
%   XYZ2Cam - given array for XYZ2Cam transform
%   wbcoeffs - a vector containing the R, G, B coefficients for white balancing.

    % read RAW image
    warning off MATLAB:tifflib:TIFFReadDirectory:libraryWarning
    obj = Tiff(filename, 'r');
    offsets = getTag(obj, 'SubIFD');
    setSubDirectory(obj, offsets(1));
    rawimg = read(obj);
    close(obj);

    % get metadata
    meta_info = imfinfo(filename);
    % (x_origin, y_origin) is the uper left corner of the useful part of
    % the sensor and consequently of the array rawim
    y_origin = meta_info.SubIFDs{1}.ActiveArea(1)+1;
    x_origin = meta_info.SubIFDs{1}.ActiveArea(2)+1;
    % width and height of the image (the useful part of array rawim)
    width = meta_info.SubIFDs{1}.DefaultCropSize(1);
    height = meta_info.SubIFDs{1}.DefaultCropSize(2);

    blacklevel = meta_info.SubIFDs{1}.BlackLevel(1); % sensor level corresponding to black
    whitelevel = meta_info.SubIFDs{1}.WhiteLevel; % sensor level corresponding to white

    wbcoeffs = (meta_info.AsShotNeutral).^-1;
    wbcoeffs = wbcoeffs/wbcoeffs(2); % green channel will be left unchanged

    XYZ2Cam = meta_info.ColorMatrix2;
    XYZ2Cam = reshape(XYZ2Cam,3,3)';

    % Crop to only valid pixels
    rawimg = double(rawimg(y_origin:y_origin+height-1, x_origin:x_origin+width-1));

    % Set black to 0 and white to 1 and crop values outside [0, 1].
    lin_bayer = (rawimg-blacklevel)/(whitelevel-blacklevel);
    lin_bayer = max(0, min(lin_bayer,1));

    rawimg = lin_bayer;
end
