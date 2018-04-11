clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

setup_par

if ~isempty(strfind(param.hostname, 'Justins-Mac'))
	param.ppath = '/Users/justin/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/fish/20180227'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/fish/20180227'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
 	param.ipath = '/om/project/boyden/DualLensLightField/01_23_18_NLS/video1';
	param.opath = '/om/user/jkinney/DLFM/fish/20180227/vid1';
end
param.inputFilePath1 = [param.ipath '/vid1/interpolate/'];
param.inputFilePath2 = [param.ipath '/confocal_data/'];

%%
param.savePath = [param.opath '/interpolate/'];
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 0.5;     % um
param.interp = 1;

param.angle   = [-1.0*pi/2 0 pi]; % radians
param.Nnull = 100;
param.confocal = true;

param

if 0>1
    param.inputFileName = {'Recon3D_solver_1_FrameNumber_0001.mat'};
    register(param)
    % close all;
    %
    % param.inputFileName = {'Recon3D_solver_1_FrameNumber_0500.mat'};
    % register(param)
    % close all;
else
    param.rapid = true;
    param.inputFileName = {'Recon3D_solver_1_FrameNumber_0001.mat'};
    % from /Users/justin/Desktop/DDLFM/fish/20180227/interpolate/
    % Recon3D_solver_1_FrameNumber_0001__20180407_131819.log
    param.trans = [302.8636 192.9936 77.8453];
    param.rot = [-1.4239 0.0069 3.1928];
    %param.centroid = [298.4520 298.4520 142];
    register(param)
end







