clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

%setup_par

if ~isempty(strfind(param.hostname, 'Justins-Mac'))
	param.ppath = '/Users/justin/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180602'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180602'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = '/om/user/jkinney/DLFM/bead/20180602';
	param.opath = param.ipath;
        %param.opath = '/om/scratch/Tue/jkinney/bead/20180602';
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

%%
param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.5; % um
param.voxel_y = 0.5; % um
param.voxel_z = 4.0;     % um

param.angle   = [pi/2 0 0]; % radians
%param.rot_amp = [0.3 0.1 0.1]; 

param.Pmelt = 0.01;
param.trans_amp = 2;
param.rot_amp = [2 2 2];
param.T0 = 1e9;
param.Tmin = 1e1;
param.Nnull = 1000;
param.T_fast = true;
param.frozen2 = 10000;

param

i=9;

param.inputFileName = {['Recon3D_after_' num2str(i) '_Mono_N15.mat']};
register(param)
close all;

i=10;

param.inputFileName = {['Recon3D_before_' num2str(i) '_Mono_N15.mat']};
register(param)
close all;

param.inputFileName = {['Recon3D_after_' num2str(i) '_Mono_N15.mat']};
register(param)
close all;

% for i=5:10
%     param.inputFileName = {['Recon3D_before_' num2str(i) '_Mono_N15.mat']};
%     register(param)
%     close all;
%
%     param.inputFileName = {['Recon3D_after_' num2str(i) '_Mono_N15.mat']};
%     register(param)
%     close all;
% end
