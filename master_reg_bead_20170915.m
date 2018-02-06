clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);



if ~isempty(strfind(param.hostname, 'Justins-Mac'))
	%% setup parallel pool
	delete(gcp('nocreate'));
	parpool(2);


	param.ppath = '/Users/justin/Desktop/DLFM'
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20170915']
	param.opath = param.ipath

elseif ~isempty(strfind(param.hostname, 'willis'))
	%% setup parallel pool
	delete(gcp('nocreate'));
	parpool(2);

	param.ppath = '/home/jkinney/Desktop/DLFM'
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20170915']
	param.opath = param.ipath
else
	%% setup parallel pool
	delete(gcp('nocreate'));
	parpool;

	param.ppath = '/om/user/jkinney/DLFM'
	addpath([param.ppath '/lotus-registration']);
    	param.ipath = '/om/user/jkinney/DLFM/bead/20170915'
	param.opath = [param.ppath '/bead/20170915']
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0;     % um
param.clip    = [100,100,100,100,0,0];
param.scale_trans = 40;
param.scale_rot   = 20;
param.trans_amp = param.scale_trans * param.voxel_x; % um
param.rot_amp = param.scale_rot * pi/800; % radians

param

param.inputFileName = {'raw.mat'};
register(param)
close all;