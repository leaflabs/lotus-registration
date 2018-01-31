clear all;
close all;

%%
[~,hostname] = system('hostname')
hostname = strtrim(hostname);

if strcmp(hostname,'Justins-Mac.local')

	ppath = '/Users/justin/Desktop/DLFM'
	addpath(ppath);
	ipath = [ppath '/worm/20171215']
	opath = ipath

elseif strcmp(hostname,'willis')

	ppath = '/home/jkinney/Desktop/DLFM'
	addpath(ppath);
	ipath = [ppath '/worm/20171215']
	opath = ipath
else
	ppath = '/om/user/jkinney/DLFM'
	addpath(ppath);
	ipath = '/om/user/ehoseini/MyData/DuallensLightField/7_20_17/video_1'
	opath = [ppath '/worm/20171215']
end
inputFilePath1 = [ipath '/horizontal/Reconstructed/'];
inputFilePath2 = [ipath '/vertical/Reconstructed/'];

%%
savePath = [opath '/interpolate/'];
voxel_x = 0.323/2; % um
voxel_y = 0.323/2; % um
voxel_z = 2.0; % um
angle   = [-1.0*pi/2 0 0];
clip    = [450,250,150,350,0,0];
myfunc_combine  = 'multiply_sqrt';
myfunc_MI  = 'multiply';
scale_trans = 40;
scale_rot   = 20;

plot = true;

%%
offset  = [-5 -5 0];
inputFileName = {'Recon3D_solver_1_FrameNumber_0001.mat'};
register

close all;
clearvars -except angle inputFilePath2  myfunc ppath scale_rot voxel_x ...
    voxel_z clip inputFilePath1 ipath opath savePath scale_trans voxel_y;

offset  = [-5 -5 0];
inputFileName = {'Recon3D_solver_1_FrameNumber_1000.mat'};
register
