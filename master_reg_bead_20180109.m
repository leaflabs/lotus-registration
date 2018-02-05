clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

if ~isempty(strfind(param.hostname, 'Justins-Mac'))

	param.ppath = '/Users/justin/Desktop/DLFM'
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180109']
	param.opath = param.ipath

elseif ~isempty(strfind(param.hostname, 'willis'))

	param.ppath = '/home/jkinney/Desktop/DLFM'
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180109']
	param.opath = param.ipath
else
	param.ppath = '/om/user/jkinney/DLFM'
	addpath([param.ppath '/lotus-registration']);
    param.ipath = '/om/user/jkinney/DLFM/bead/20180109'
	param.opath = [param.ppath '/worm/20171215']
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323/2; % um
param.voxel_y = 0.323/2; % um
param.voxel_z = 2.0;     % um
param.clip    = [100,100,100,100,0,0];
param.scale_trans = 40;
param.scale_rot   = 20;
param.trans_amp = param.scale_trans * param.voxel_x; % um
param.rot_amp = param.scale_rot * pi/800; % radians

param.inputFileName = {'Recon3D_1_Mono_N15.mat'};

param

register(param)
close all;

param.inputFileName = {'Recon3D_2_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_3_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_4_Mono_N15.mat'};
register(param)
close all;

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

