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
	param.ipath = [param.ppath '/bead/20180328'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180328'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = '/om/user/jkinney/DLFM/bead/20180328';
	param.opath = param.ipath;
        %param.opath = '/om/scratch/Tue/jkinney/bead/20180328';
    param.parallel = false;
end
%param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath1 = [param.ipath '/vertical_phantom/'];

%param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical_phantom/copies_of_vertical/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.5; % um
param.voxel_y = 0.5; % um
param.voxel_z = 4.0;     % um

param.angle   = [0 0 0]; % radians
param.interp = 1;

param

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_0_0_trans_0_0_50.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_0_0_trans_0_50_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_0_0_trans_50_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_0_0_trans_50_50_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_0_0_trans_50_50_50.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_0_90_trans_0_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_0_90_0_trans_0_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_90_0_0_trans_0_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_90_0_90_trans_0_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_90_90_0_trans_0_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_-90_-90_-90_trans_0_0_0.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_after_9_Mono_N15_20180608_110054_rot_90_90_90_trans_0_0_0.mat'};
register(param)
close all;

