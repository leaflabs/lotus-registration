clear all;
close all;

%%
[~,hostname] = system('hostname')
hostname = strtrim(hostname);

if strcmp(hostname,'Justins-Mac.local')

	ppath = '/Users/justin/Desktop/DLFM'
	addpath([ppath '/lotus-registration']);
	ipath = [ppath '/bead/20171215']
	opath = ipath

elseif strcmp(hostname,'willis')

	ppath = '/home/jkinney/Desktop/DLFM'
	addpath([ppath '/lotus-registration']);
	ipath = [ppath '/bead/20171215']
	opath = ipath
else
	ppath = '/om/user/jkinney/DLFM'
	addpath([ppath '/lotus-registration']);
    ipath = '/om/user/jkinney/DLFM/bead/20180109'
	opath = ipath
end
inputFilePath1 = [ipath '/horizontal/Reconstructed/'];
inputFilePath2 = [ipath '/vertical/Reconstructed/'];

%%
savePath = [opath '/interpolate/'];
voxel_x = 0.323/2; % um
voxel_y = 0.323/2; % um
voxel_z = 2.0; % um
clip    = [100,100,100,100,0,0];

%
scale_trans = 40;
scale_rot   = 20;
inputFileName = {'Recon3D_1_Mono_N15.mat'};
register

close all;
clearvars -except ppath ipath opath...
   inputFilePath1 inputFilePath2 savePath...
   voxel_x voxel_y voxel_z...
   angle clip scale_trans...
   myfunc_combine myfunc_MI...
   ;

scale_trans = 4;
scale_rot   = 2;
inputFileName = {'Recon3D_solver_1_FrameNumber_0001.mat'};
register

 close all;
 clearvars -except ppath ipath opath...
     inputFilePath1 inputFilePath2 savePath...
     voxel_x voxel_y voxel_z...
     angle clip scale_trans...
     myfunc_combine myfunc_MI...
     ;

scale_trans = 40;
scale_rot   = 20;
inputFileName = {'Recon3D_solver_1_FrameNumber_1000.mat'};
register

close all;
clearvars -except ppath ipath opath...
   inputFilePath1 inputFilePath2 savePath...
   voxel_x voxel_y voxel_z...
   angle clip scale_trans...
   myfunc_combine myfunc_MI...
   ;

scale_trans = 4;
scale_rot   = 2;
inputFileName = {'Recon3D_solver_1_FrameNumber_1000.mat'};
register
