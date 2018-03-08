clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

setup_par

if ~isempty(strfind(param.hostname, 'Justins-Mac'))
	param.ppath = '/Users/justin/Desktop/DLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180109'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180109'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
    param.ipath = '/om/user/jkinney/DLFM/bead/20180109';
	param.opath = param.ipath;
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323/2; % um
param.voxel_y = 0.323/2; % um
param.voxel_z = 2.0;     % um

param

param.inputFileName = {'Recon3D_1_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_2_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_3_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_4_Mono_N15.mat'};
tmp = param.angle;
%param.angle   = [-1.2*pi/2 0 0];
register(param)
close all;
param.angle   = tmp;

param.inputFileName = {'Recon3D_5_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_6_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_7_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_8_Mono_N15.mat'};
register(param)
close all;

