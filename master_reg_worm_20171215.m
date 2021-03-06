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
	param.ipath = [param.ppath '/worm/20171215'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/worm/20171215'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
 	param.ipath = '/om/project/boyden/DualLensLightField/12_15_17_40X_worm/video 1';
	param.opath = '/om/user/jkinney/DLFM/worm/20171215';
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/interpolate/'];
param.voxel_x = 0.323/2; % um
param.voxel_y = 0.323/2; % um
param.voxel_z = 2.0;     % um

%param.angle   = [-0.5*pi/2 0 0];

param

param.inputFileName = {'Recon3D_solver_1_FrameNumber_0001.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_solver_1_FrameNumber_1000.mat'};
register(param)
close all;

