clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);



if ~isempty(strfind(param.hostname, 'Justins-Mac'))
	%% setup parallel pool
	p = gcp('nocreate');
    	if isempty(p)
        	parpool(2);
    	elseif ~p.NumWorkers==2
        	delete(p);
        	parpool(2);
        	p = gcp('nocreate');
    	end


	param.ppath = '/Users/justin/Desktop/DLFM'
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180109']
	param.opath = param.ipath

elseif ~isempty(strfind(param.hostname, 'willis'))
	%% setup parallel pool
	p = gcp('nocreate');
    	if isempty(p)
        	parpool(2);
    	elseif ~p.NumWorkers==2
        	delete(p);
        	parpool(2);
        	p = gcp('nocreate');
    	end

	param.ppath = '/home/jkinney/Desktop/DLFM'
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180109']
	param.opath = param.ipath
else
	%% setup parallel pool
	p = gcp;

	param.ppath = '/om/user/jkinney/DLFM'
	addpath([param.ppath '/lotus-registration']);
    	param.ipath = '/om/user/jkinney/DLFM/bead/20180109'
	param.opath = [param.ppath '/bead/20180109']
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
param.trans_amp = param.scale_trans * param.voxel_x; % um
param.rot_amp = param.scale_rot * pi/800; % radians

param

param.inputFileName = {'Recon3D_1_Mono_N15.mat'};
%param.offset = {[0], [], []};
register(param)
close all;
%param.offset = {[], [], []};

param.inputFileName = {'Recon3D_2_Mono_N15.mat'};
%param.offset = {[], [-30], []};
register(param)
close all;
%param.offset = {[], [], []};

param.inputFileName = {'Recon3D_3_Mono_N15.mat'};
register(param)
close all;

param.inputFileName = {'Recon3D_4_Mono_N15.mat'};
param.angle   = [-1.2*pi/2 0 0];
register(param)
close all;
param.angle   = [0 0 0];

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

