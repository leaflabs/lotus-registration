clear all;
close all;

%%
[~,hostname] = system('hostname')
hostname = strtrim(hostname);

if strcmp(hostname,'Justins-Mac.local')

	ppath = '/Users/justin/Desktop/DLFM'
	addpath([ppath '/lotus-registration']);
	ipath = [ppath '/worm/20171215']
	opath = ipath

elseif strcmp(hostname,'willis')

	ppath = '/home/jkinney/Desktop/DLFM'
	addpath([ppath '/lotus-registration']);
	ipath = [ppath '/beads/20180109']
	opath = ipath
else
	ppath = '/om/user/jkinney/DLFM'
	addpath([ppath '/lotus-registration']);
    ipath = '/om/project/boyden/DualLensLightField/12_15_17_40X_worm/video 1'
	opath = [ppath '/worm/20171215']
end
inputFilePath1 = [ipath '/horizontal/Reconstructed/'];
inputFilePath2 = [ipath '/vertical/Reconstructed/'];

%% Recon3D_1
% FAILED - LFM1 and LFM2 not same scale, maybe
savePath = [opath '/registration/'];
voxel_x = 0.323/2; % um
voxel_y = 0.323/2; % um
voxel_z = 2.0; % um
clip    = [100,500,100,100,0,0];
scale_trans = 40;
scale_rot   = 20;
offset = {[], [], []};
inputFileName = {'Recon3D_1_Mono_N15.mat'};
register

%% Recon3D_2

% savePath = [opath '/registration/'];
% voxel_x = 0.323/2; % um
% voxel_y = 0.323/2; % um
% voxel_z = 2.0; % um
% clip    = [100,100,100,100,0,0];
% scale_trans = 4;
% scale_rot   = 2;
% offset = {[], [-40], []};
% inputFileName = {'Recon3D_2_Mono_N15.mat'};
% register
% 
% close all;
% clearvars -except ppath ipath opath...
%    inputFilePath1 inputFilePath2 savePath...
%    voxel_x voxel_y voxel_z...
%    angle clip scale_trans...
%    myfunc_combine myfunc_MI...
%    ;
% 
% savePath = [opath '/registration/'];
% voxel_x = 0.323/2; % um
% voxel_y = 0.323/2; % um
% voxel_z = 2.0; % um
% clip    = [100,100,100,100,0,0];
% scale_trans = 40;
% scale_rot   = 20;
% offset = {[], [], []};
% inputFileName = {'Recon3D_2_Mono_N15.mat'};
% register
% 
% close all;
% clearvars -except ppath ipath opath...
%    inputFilePath1 inputFilePath2 savePath...
%    voxel_x voxel_y voxel_z...
%    angle clip scale_trans...
%    myfunc_combine myfunc_MI...
%    ;
% 
% savePath = [opath '/registration/'];
% voxel_x = 0.323/2; % um
% voxel_y = 0.323/2; % um
% voxel_z = 2.0; % um
% clip    = [100,100,100,100,0,0];
% scale_trans = 40;
% scale_rot   = 20;
% offset = {[], [-40], []};
% inputFileName = {'Recon3D_2_Mono_N15.mat'};
% register
% 
% close all;
% clearvars -except ppath ipath opath...
%    inputFilePath1 inputFilePath2 savePath...
%    voxel_x voxel_y voxel_z...
%    angle clip scale_trans...
%    myfunc_combine myfunc_MI...
%    ;





