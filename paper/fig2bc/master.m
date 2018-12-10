close all;
clear all;

config = struct;
config.fpath = '/Users/justin/Desktop/DDLFM/fig2bc/mat_files/';
config.xypath = config.fpath;
config.outpath = '/Users/justin/Desktop/DDLFM/fig2bc/out/';
config.tmpf = [config.outpath 'bead_summary_data.mat'];

config.timestamp = datestr(datetime('now'),'yyyymmdd_HHMMSS');
fdiary = sprintf('%s%s.log',config.outpath,config.timestamp);
diary(fdiary)
disp('%%');
disp(['%% ' fdiary]);
disp('%%');

config.div = 4;
config.pixel = 0.5; % um
config.figpos = [1 1 1662 976];

if ~exist(config.tmpf,'file')
    
    bead_struct.header = ['header = index '...
        'centroid_dim1(um) centroid_dim2(um) centroid_dim3(um) '...
        'fwhm_dim1a(um) fwhm_dim1b(um) fwhm_dim2(um) fwhm_dim3(um)'];
    
    config.fname = '04-Oct-2017 17_07_02_multiply.mat';
    config.xyname = 'xycoords_side2.mat';
    config.zspacing = 0.5;
    config.crop_um = 20;
    config.crop_z = ceil(config.crop_um/config.zspacing);
    config.crop_row = ceil(config.crop_um/config.pixel);
    config.label = 'DLFM_';
    config
    bead_psf = axial_psf(config);
    bead_struct.dlfm = bead_psf;
    
    config.fname = 'raw_side2.mat';
    config.xyname = 'xycoords_side2.mat';
    config.zspacing = 4;
    config.crop_um = 50;
    config.crop_z = ceil(config.crop_um/config.zspacing);
    config.crop_row = ceil(config.crop_um/config.pixel);
    config.label = 'LFM2_';
    config
    bead_psf = axial_psf(config);
    bead_struct.lfm2 = bead_psf;
    
    config.fname = 'raw_side1.mat';
    config.xyname = 'xycoords_side1.mat';
    config.zspacing = 4;
    config.crop_um = 50;
    config.crop_z = ceil(config.crop_um/config.zspacing);
    config.crop_row = ceil(config.crop_um/config.pixel);
    config.label = 'LFM1_';
    config
    bead_psf = axial_psf(config);
    bead_struct.lfm1 = bead_psf;
    
    save(config.tmpf,'config','bead_struct');
else
    load(config.tmpf);
end


%%% lists

% black
a = find(bead_struct.dlfm(:,2)<0)';
config.blacklist_dlfm =  sort([a [14 16 18 19 34 35 43 44 45 53]]);
a = find(bead_struct.lfm2(:,2)<0)';
config.blacklist_lfm2 =  sort([a [12 14 16 18 19 24 25 31 32 35 38 39 44 46 47 50 52 53 54 56 61 66]]);
a = find(bead_struct.lfm1(:,2)<0)';
config.blacklist_lfm1 =  sort([a [2 5 6 14 37 42 45 46 47 48 49]]);

% white
a = size(bead_struct.dlfm);
config.whitelist_dlfm = setdiff(1:a(1),config.blacklist_dlfm);
a = size(bead_struct.lfm2);
config.whitelist_lfm2 = setdiff(1:a(1),config.blacklist_lfm2);
a = size(bead_struct.lfm1);
config.whitelist_lfm1 = setdiff(1:a(1),config.blacklist_lfm1);

%%% sanity checks

% todo: to fix this code, identify overlap in whitelist for dlfm and lfm2
%
% centroid of bead in two volumes should be similar
%centroid_diff = [bead_struct.dlfm(:,1) bead_cell{2}(:,2:4)- bead_cell{3}(:,2:4)];

% for each bead fwhm1 should be similar when measured twice - fwhm1a and fwhm1b
wl = config.whitelist_dlfm;
fwhm1_diff_dlfm = [bead_struct.dlfm(wl,1) bead_struct.dlfm(wl,5)-bead_struct.dlfm(wl,6)]
wl = config.whitelist_lfm2;
fwhm1_diff_lfm2 = [bead_struct.lfm2(wl,1) bead_struct.lfm2(wl,5)-bead_struct.lfm2(wl,6)]
wl = config.whitelist_lfm1;
fwhm1_diff_lfm1 = [bead_struct.lfm1(wl,1) bead_struct.lfm1(wl,5)-bead_struct.lfm1(wl,6)]

%%% resolution

% compare the mean half width across the conditions (LFM1, LFM2, DLF), for y and z.
fprintf('\nhalf width (um)\n[dim1 dim2 dim3] \n');
print_ms( bead_struct.lfm2, config.whitelist_lfm2, 'LFM2');
print_ms( bead_struct.lfm1, config.whitelist_lfm1, 'LFM1');
print_ms( bead_struct.dlfm, config.whitelist_dlfm, 'DLFM');

% todo: model a pair of neaighboring beads as two gaussians
%       in order to calculate their separability
%       as a function of distance between the beads.

%%% isotropy

% compare the ratio of y:z across the conditions.
fprintf('\nratio of fwhm dim1:dim3 and dim2:dim3\n');
print_fwhm_ratio( bead_struct.lfm2, config.whitelist_lfm2, 'LFM2');
print_fwhm_ratio( bead_struct.lfm1, config.whitelist_lfm1, 'LFM1');
print_fwhm_ratio( bead_struct.dlfm, config.whitelist_dlfm, 'DLFM');

% todo: look at the distribution of the ratios, ie is it normaly distributed,
%       and is the mean different from 1 (based on null hypothesis you will be
%       able to say something about how isotropic the resolution is). 

diary off

