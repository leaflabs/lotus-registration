close all;
clear all;

config = struct;
fvec = {'04-Oct-2017 17_07_02_multiply.mat'};
config.fpath = '/Users/justin/Desktop/DDLFM/fig2bc/mat_files/';
config.xypath = config.fpath;
config.outpath = '/Users/justin/Desktop/DDLFM/fig2bc/out/';
config.tmpf = [config.outpath 'bead_fwhm.mat'];

timestamp = datestr(datetime('now'));
fdiary = sprintf('%s%s.log',config.outpath,timestamp);
diary(fdiary)
disp('%%');
disp(['%% ' datestr(datetime)]);
disp('%%');

config.div = 4;
config.pixel = 0.5; % um
config.figpos = [1           1        1662         976];

config.blacklist = [6 8 11 12 14 15 16 18 19 24 25 27 31 32 34:40 43:57 61 62 64:67];

if ~exist(config.tmpf,'file')
    
    bead_cell = {'header = index centroid_dim1(um) centroid_dim2(um) centroid_dim3(um) fwhm_dim1a(um) fwhm_dim1b(um) fwhm_dim2(um) fwhm_dim3(um)'};
    
    config.fname = fvec{1};
    config.zspacing = 0.5;
    config.xyname = 'xycoords_side2.mat';
    config.crop_um = 20;
    config.crop_z = ceil(config.crop_um/config.zspacing);
    config.crop_row = ceil(config.crop_um/config.pixel);
    config
    bead_psf = axial_psf(config);
    bead_cell = [bead_cell {bead_psf}];
    
    config.fname = 'raw_side2.mat';
    config.zspacing = 4;
    config.xyname = 'xycoords_side2.mat';
    config.crop_um = 50;
    config.crop_z = ceil(config.crop_um/config.zspacing);
    config.crop_row = ceil(config.crop_um/config.pixel);
    config
    bead_psf = axial_psf(config);
    bead_cell = [bead_cell {bead_psf}];
    
    % fname = 'raw_side1.mat';
    % xyname = 'xycoords_side1.mat';
    % close all;
    % axial_psf
    % bead_cell = [bead_cell {bead_psf}];
    
    save(config.tmpf,'config','bead_cell');
else
    load(config.tmpf);
end

%%% sanity checks
% centroid of bead in two volumes should be similar
centroid_diff = [bead_cell{2}(:,1) bead_cell{2}(:,2:4)- bead_cell{3}(:,2:4)];
% for each bead fwhm1 should be similar when measured twice - fwhm1a and fwhm1b
fwhm1_diff = [bead_cell{2}(:,1) bead_cell{2}(:,5)- bead_cell{3}(:,6) ...
    bead_cell{3}(:,5)- bead_cell{3}(:,6)];

% compare the mean half width across the conditions (LFM1, LFM2, DLF), for y and z.
fprintf('\nhalf width (um)\n[dim1 dim2 dim3] \n');
tmp = mean(bead_cell{3}(:,6:8));
fprintf('LFM2 mean [%2.1f %2.1f %2.1f]\n',tmp(1),tmp(2),tmp(3));
tmp = std(bead_cell{3}(:,6:8));
fprintf('     stddev [%2.1f %2.1f %2.1f]\n',tmp(1),tmp(2),tmp(3));
tmp = mean(bead_cell{2}(:,6:8));
fprintf('DLFM mean [%2.1f %2.1f %2.1f]\n',tmp(1),tmp(2),tmp(3));
tmp = std(bead_cell{2}(:,6:8));
fprintf('     stddev [%2.1f %2.1f %2.1f]\n',tmp(1),tmp(2),tmp(3));

% compare the ratio of y:z across the conditions.
fprintf('\nratio of fwhm dim1:dim3 and dim2:dim3\n');
tmp1 = mean(bead_cell{3}(:,6)./bead_cell{3}(:,8));
tmp2 = mean(bead_cell{3}(:,7)./bead_cell{3}(:,8));
fprintf('LFM2 mean [1:3__%2.1f, 2:3__%2.1f]\n',tmp1,tmp2);
tmp1 = std(bead_cell{3}(:,6)./bead_cell{3}(:,8));
tmp2 = std(bead_cell{3}(:,7)./bead_cell{3}(:,8));
fprintf('     stddev [1:3__%2.1f, 2:3__%2.1f]\n',tmp1,tmp2);
tmp1 = mean(bead_cell{2}(:,6)./bead_cell{2}(:,8));
tmp2 = mean(bead_cell{2}(:,7)./bead_cell{2}(:,8));
fprintf('DLFM mean [1:3__%2.1f, 2:3__%2.1f]\n',tmp1,tmp2);
tmp1 = std(bead_cell{2}(:,6)./bead_cell{2}(:,8));
tmp2 = std(bead_cell{2}(:,7)./bead_cell{2}(:,8));
fprintf('     std [1:3__%2.1f, 2:3__%2.1f]\n',tmp1,tmp2);

% could also look at the distribution of the ratios, ie is it normaly distributed,
% and is the mean different from 1 (based on null hypothesis you will be
% able to say something about how isotropic the resolution is). 

diary off

