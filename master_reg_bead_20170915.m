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
	param.ipath = [param.ppath '/bead/20170915'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20170915'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
    param.ipath = '/om/user/jkinney/DLFM/bead/20170915';
	param.opath = param.ipath;
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0;     % um

%param.angle   = [-1.2*pi/2 0 0];

param.inputFileName = {'raw.mat'};
register(param)
close all;
