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
	param.ipath = [param.ppath '/bead/20180406'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180406'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = '/om/user/jkinney/DLFM/bead/20180406';
	param.opath = param.ipath;
        %param.opath = '/om/scratch/Tue/jkinney/bead/20180406';
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0;     % um

param.angle   = [-1.2*pi/2 0 0]; % radians
param.rot_amp = [0.3 0.1 0.1]; 

param

param.inputFileName = {'Recon3D_camera_ccw_1_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_camera_ccw_2_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_camera_cw_1_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_camera_cw_2_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_minus_1_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_minus_2_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_minus_3_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_minus_4_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_plus_1_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_plus_2_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_plus_3_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_micro_plus_4_Mono_N15.mat'};
register(param)
close all;
param.inputFileName = {'Recon3D_standard_Mono_N15.mat'};
register(param)
close all;
