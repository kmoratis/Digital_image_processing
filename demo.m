% Demo script to display the readdng() and dng2rgb() functions

% using readdng() to read the RAW image and the useful metadata
[rawim, XYZ2Cam, wbcoeffs] = readdng('./RawImage.DNG');

% using dng2rgb to transform RAW image (.DNG) to RGB
M = 4000;
N = 6000;
bayertype = 'RGGB';
method = 'linear';
[Csrgb,Clinear,Cxyz,Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, bayertype, method, M, N);

% create a histogram for every channel of every picture
% Ccam
R = Ccam(:,:,1);
G = Ccam(:,:,2);
B = Ccam(:,:,3);

f = figure('Name', 'Ccam Histograms');
f.Position(3:4) = [1200 600]; % Widht, height (px)
f.Position(1:2) = [100 100];
t = tiledlayout(1,3,'TileSpacing','Compact');

% R channel
nexttile
histogram(R, 'Normalization','probability','NumBins',100,'FaceColor','red');
% G channel
nexttile
histogram(G, 'Normalization','probability','NumBins',100,'FaceColor','green');
% B channel
nexttile
histogram(B, 'Normalization','probability','NumBins',100,'FaceColor','blue');

title(t, 'Histograms for Image 1');
ylabel(t, 'Relative frequency');

% Cxyz
R = Cxyz(:,:,1);
G = Cxyz(:,:,2);
B = Cxyz(:,:,3);

f = figure('Name', 'Cxyz Histograms');
f.Position(3:4) = [1200 600]; % Widht, height (px)
f.Position(1:2) = [100 100];
t = tiledlayout(1,3,'TileSpacing','Compact');

% R channel
nexttile
histogram(R, 'Normalization','probability','NumBins',100,'FaceColor','red');
% G channel
nexttile
histogram(G, 'Normalization','probability','NumBins',100,'FaceColor','green');
% B channel
nexttile
histogram(B, 'Normalization','probability','NumBins',100,'FaceColor','blue');

title(t, 'Histograms for Image 2');
ylabel(t, 'Relative frequency');

% Clinear
R = Clinear(:,:,1);
G = Clinear(:,:,2);
B = Clinear(:,:,3);

f = figure('Name', 'Clinear Histograms');
f.Position(3:4) = [1200 600]; % Widht, height (px)
f.Position(1:2) = [100 100];
t = tiledlayout(1,3,'TileSpacing','Compact');

% R channel
nexttile
histogram(R, 'Normalization','probability','NumBins',100,'FaceColor','red');
% G channel
nexttile
histogram(G, 'Normalization','probability','NumBins',100,'FaceColor','green');
% B channel
nexttile
histogram(B, 'Normalization','probability','NumBins',100,'FaceColor','blue');

title(t, 'Histograms for Image 3');
ylabel(t, 'Relative frequency');

% Csrgb
R = Csrgb(:,:,1);
G = Csrgb(:,:,2);
B = Csrgb(:,:,3);

f = figure('Name', 'Csrgb Histograms');
f.Position(3:4) = [1200 600]; % Widht, height (px)
f.Position(1:2) = [100 100];
t = tiledlayout(1,3,'TileSpacing','Compact');

% R channel
nexttile
histogram(R, 'Normalization','probability','NumBins',100,'FaceColor','red');
% G channel
nexttile
histogram(G, 'Normalization','probability','NumBins',100,'FaceColor','green');
% B channel
nexttile
histogram(B, 'Normalization','probability','NumBins',100,'FaceColor','blue');

title(t, 'Histograms for Image 4');
ylabel(t, 'Relative frequency');
