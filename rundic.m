addpath(genpath('FIDIC'));
sSize = [64 64];
incORcum = 'c'; %use 'i' for incremental mode and 'c' for cumulative
norm_xcc = 'u'; %use 'norm' for normalized cross-correlation, considerable time-cost
ext_in = '*.jpg'; %Input image format

max_def_idx = 'b'; %Specify where the max deformation occurs
yn = 'y';

images=cell(2,1);
images{1}=I;
images{2}=I1;
tic
[u, cc, dm] = funIDIC(images, sSize, incORcum, norm_xcc);
toc