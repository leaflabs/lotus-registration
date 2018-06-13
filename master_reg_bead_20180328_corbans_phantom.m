clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

%setup_par

if ~isempty(strfind(param.hostname, 'Justins-Mac')) || ~isempty(strfind(param.hostname, 'dhcp-18-189-20-128.dyn.mit.edu'))
	param.ppath = '/Users/justin/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = param.ppath;
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = param.ppath;
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = '/om/user/jkinney/DLFM';
	param.opath = param.ipath;
    param.parallel = false;
end
%param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath1 = [param.ipath '/Transformation-Phantoms-master/analysis/'];

%param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];
param.inputFilePath2 = [param.ipath ...
    '/Transformation-Phantoms-master/analysis/copies_of_rot_0_0_0_trans_0_0_0/'];

%%
param.savePath = [param.inputFilePath1 '/registration/'];
param.voxel_x = 1; % um
param.voxel_y = 1; % um
param.voxel_z = 1;     % um

param.angle   = [0 0 0]; % radians
param.interp = 1;
param.clip = [0 0 0 0 0 0];
param.confocal = true; % to rotate in all 3 dims.
param.Pmelt = 0.005;
param.trans_amp = 1;
param.rot_amp = [1 1 1];

% %% new set

param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_0_0_trans_0_0_0.mat'};

% param.angle   = [0*pi/180 0 0]; % radians
% register(param)
% close all;
% 
% param.angle   = [15*pi/180 0 0]; % radians
% register(param)
% close all;
% 
% param.angle   = [30*pi/180 0 0]; % radians
% register(param)
% close all;
% 
% param.angle   = [45*pi/180 0 0]; % radians
% register(param)
% close all;
% 
% param.angle   = [60*pi/180 0 0]; % radians
% register(param)
% close all;
% 
% param.angle   = [75*pi/180 0 0]; % radians
% register(param)
% close all;

% param.angle   = [90*pi/180 0 0]; % radians
% register(param)
% close all;
% 
% param.angle   = [105*pi/180 0 0]; % radians
% register(param)
% close all;

% param.angle   = [120*pi/180 0 0]; % radians
% register(param)
% close all;

% param.angle   = [135*pi/180 0 0]; % radians
% register(param)
% close all;

%% original set

angle_0 = 85; % degrees

param.frozen2 = 1000;

% param.angle   = [angle_0/180*pi angle_0/180*pi 0]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_90_0_trans_0_0_0.mat'};
% register(param)
% close all;
% 
% param.angle   = [0 angle_0/180*pi 0]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_90_0_trans_0_0_0.mat'};
% register(param)
% close all;
% 
% param.angle   = [0 0 angle_0/180*pi]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_90_trans_0_0_0.mat'};
% register(param)
% close all;
% 
% param.angle   = [angle_0/180*pi 0 0]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_0_0_trans_0_0_0.mat'};
% register(param)
% close all;
% 
% param.angle   = [angle_0/180*pi 0 angle_0/180*pi]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_0_90_trans_0_0_0.mat'};
% register(param)
% close all;

param.angle   = [-angle_0/180*pi -angle_0/180*pi -angle_0/180*pi]; % radians
param.inputFileName = {'axis_phantom_100_20180612_201406_rot_-90_-90_-90_trans_0_0_0.mat'};
register(param)
close all;

param.angle   = [angle_0/180*pi angle_0/180*pi angle_0/180*pi]; % radians
param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_90_90_trans_0_0_0.mat'};
register(param)
close all;

keyboard

param.angle   = [0 0 0]; % radians

param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_0_0_50.mat'};
register(param)
close all;

param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_0_50_0.mat'};
register(param)
close all;

param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_50_0_0.mat'};
register(param)
close all;

param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_50_50_0.mat'};
register(param)
close all;

param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_50_50_50.mat'};
register(param)
close all;



