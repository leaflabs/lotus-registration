clear all;
close all;

%% load parameters
param = init_param();

%% setup parameters
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

setup_par

if ~isempty(strfind(param.hostname, 'Justins-Mac'))
        param.ppath = '/Users/justin/Desktop/DLFM';
        addpath([param.ppath '/lotus-registration']);
        param.ipath = [param.ppath '/fish/20180227/vid1'];
        param.opath = param.ipath;;
        param.inter = [param.ipath '/interpolate/'];
elseif ~isempty(strfind(param.hostname, 'willis'))
        param.ppath = '/home/jkinney/Desktop/DLFM';
        addpath([param.ppath '/lotus-registration']);
        param.ipath = [param.ppath '/fish/20180227/vid1'];
        param.opath = param.ipath;
        param.inter = [param.ipath '/interpolate/'];
else
        param.ppath = '/home/jkinney';
        addpath([param.ppath '/lotus-registration']);
 	param.ipath = '/om/project/boyden/DualLensLightField/01_23_18_NLS/video1';
        param.opath = '/om/scratch/Tue/jkinney/fish/20180227/vid1';
        param.spath = '/om/user/jkinney/DLFM/fish/20180227/vid1';
        param.inter = [param.spath '/interpolate/'];

	C = strsplit(param.opath,'/');
	mypath = '';
	for i=2:length(C)
	    mypath = [mypath '/' C{i}]
	    if ~exist(mypath,'dir')
	        status = mkdir(mypath);
	        if status == 1
	            disp(['Created folder: ' mypath]);
	        else
	            disp(['Error attempting to create folder:' mypath]);
	            status
	            exit;
	        end
	    end
	end
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0;     % um

param.m = 1;
param.n = 500;
next = param.m;

param.rapid = true;

param = read_logs(param);

for i=[next:param.n]
    param.rot = param.rotM + (param.rotN-param.rotM) * (i-param.m)/(param.n-param.m);
    tmp_trans = param.transM + (param.transN-param.transM) * (i-param.m)/(param.n-param.m);
    param.trans = tmp_trans + [ param.clip(1)*param.voxel_y  ...
        param.clip(3)*param.voxel_x   param.clip(5)*param.voxel_z ];
    param.inputFileName = {sprintf('Recon3D_solver_1_FrameNumber_%04d.mat',i)};
    
    param
    
    register(param)
end


