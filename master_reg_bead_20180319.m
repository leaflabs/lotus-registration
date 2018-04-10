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
	param.ipath = [param.ppath '/bead/20180319'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180319'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = '/om/user/jkinney/DLFM/bead/20180319';
	%param.opath = param.ipath;
        param.opath = '/om/scratch/Tue/jkinney/bead/20180319';
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0;     % um

param

start = 5;

for i=start:99
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
	register(param)
	close all;
end
