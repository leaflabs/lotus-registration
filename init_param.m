function param = init_param ()
param = struct;
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0; % um
param.interp = 8;
param.myfunc_combine = 'multiply_sqrt';
param.myfunc_MI = 'multiply';

param.prate = -1;
param.Trate = 1e-1;
param.final_p = 1e-2;
param.init_p = 0.3;
param.scale_trans = 40;
param.scale_rot = [1 1 1];
param.trans_amp = param.scale_trans * param.voxel_x; % um
param.rot_amp = param.scale_rot * pi/800; % radians
%param.trans_amp = 1.0; % um
%param.rot_amp = pi/80; % radians

% amount to clip from periphery in pixels
% dim1 min, dim1 max, dim2 min, dim2 max, dim3 min, dim3 max,
param.clip = [0,0,0,0,0,0];

param.T0 = -1;
param.TC0 = 30;
param.MC0 = 30;
%param.TC0 = 1;
%param.MC0 = 900;
param.Nnull = 2000;

param.psf   = [-1.0 -1.0 -1.0]; %um
%param.offset = [-8 -30 0];
%param.offset = [0 0 0];
param.offset = {[], [], []};
param.trans = [0 0 0]; % um
%param.angle   = [-1.2*pi/2 0 0]; % radians
param.angle   = [-1.0*pi/2 0 0]; % radians
param.rot   = [0 0 0]; % radians

param.lfdisplay = false;

%param.centroid = [169.8834 184.3252 129.1565];
param.centroid = [-1 -1 -1];

% indices of side2 that are being tracked
%param.threshold = 20;
param.threshold1 = 40;
param.threshold2 = 40;
param.index2 = [];

param.N = 10000;
param.xlim_thresh = 0.9999; %HARDCODED
%param.contour_thresh = 0.99;
param.pop_thresh = 0.99; %HARDCODED
param.dynamic_range_thresh = 10;

param.rapid = false;
param.savevol = true;
param.plot = true;
%param.savePath      = '/Users/justin/Desktop/LFM volume registration/';
param.savePath      = '/home/jkinney/Desktop/LFM volume registration/registration/param_sweep/';
%param.savePath      = '/tmp/';
param.inputFilePath1 = '';
param.inputFilePath2 = '';
param.inputFileName = {
    'raw_side1.mat',
    'raw_side2.mat'};

param.xlim_1d = [0 0 0];
end

