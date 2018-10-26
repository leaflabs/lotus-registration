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
param.Pmelt = 0.01;
param.trans_amp = 2;
param.rot_amp = [2 2 2];
param.T0 = 1e6;
param.Tmin = 1e1;
% %% new set

%param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_0_0_trans_0_0_0.mat'};

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

% angle_0 = 85; % degrees
% 
% param.frozen2 = 1000;

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

% param.angle   = [-angle_0/180*pi -angle_0/180*pi -angle_0/180*pi]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_-90_-90_-90_trans_0_0_0.mat'};
% register(param)
% close all;
% 
% param.angle   = [angle_0/180*pi angle_0/180*pi angle_0/180*pi]; % radians
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_90_90_90_trans_0_0_0.mat'};
% register(param)
% close all;
% 
% param.angle   = [0 0 0]; % radians
% 
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_0_0_50.mat'};
% register(param)
% close all;
% 
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_0_50_0.mat'};
% register(param)
% close all;
% 
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_50_0_0.mat'};
% register(param)
% close all;
% 
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_50_50_0.mat'};
% register(param)
% close all;
% 
% param.inputFileName = {'axis_phantom_100_20180612_201406_rot_0_0_0_trans_50_50_50.mat'};
% register(param)
% close all;


%%
if 0>1
    angle_0 = 88; % degrees
    
    param.frozen2 = 10000;
    
    % param.trans = [75.118012 74.985333 74.737237];
    % param.rot = [1.554521 1.566939 -0.014491];
    %
    % param.trans = [74.000000 75.000000 75.000000];
    % param.rot = [1.535890 1.535890 0.000000];
    %
    % param.trans = [74.000000 75.000000 74.852057];
    % param.rot = [1.535890 1.535890 0.000000];
    %
    % param.trans = [75.512620 75.123032 75.037646];
    % param.rot = [1.556528 1.577588 -0.008215];
    %
    % param.trans = [74.877351 75.123032 75.186412];
    % param.rot = [1.555324 1.577588 -0.008215];
    %
    % param.trans = [74.877351 75.123032 74.852057];
    % param.rot = [1.539280 1.577588 -0.025145];
    %
    % param.trans = [75.832299 75.000000 74.852057];
    % param.rot = [1.539280 1.572580 -0.025145];
    %
    % param.trans = [74.000000 75.000000 74.852057];
    % param.rot = [1.535890 1.564438 0.000000];
    %
    % param.trans = [74.000000 75.000000 74.852057];
    % param.rot = [1.535890 1.564438 0.000000];
    
    % Good MI =    24772608
    % param.trans = [74.000000 75.000000 74.852057];
    % param.rot = [1.535890 1.535890 0.000000];
    
    
    % param.just_MI = false;
    % param.angle   = [angle_0/180*pi angle_0/180*pi 0]; % radians
    %     param.inputFileName = {'bead_phantom_150_20180614_083753_rot_90_90_0_trans_0_0_0.mat'};
    %     register(param)
    %     close all;
    
    %param.savePath = [param.inputFilePath1 '/registration/repeat_90_90_0_A/']; % pre-fixed bug related to last_MI
    %param.savePath = [param.inputFilePath1 '/registration/repeat_90_90_0_B/']; % fixed bug related to last_MI
    %param.savePath = [param.inputFilePath1 '/registration/repeat_90_90_0_C/']; % randomized sequence of 6 parameters
    param.savePath = [param.inputFilePath1 '/registration/repeat_90_90_0_D/']; % C + no delmi==0 allowed
    for j=1:10
        rng(j);
        j
        param.angle   = [angle_0/180*pi angle_0/180*pi 0]; % radians
        param.inputFileName = {'bead_phantom_150_20180614_083753_rot_90_90_0_trans_0_0_0.mat'};
        register(param)
        close all;
    end
    
    keyboard
    
    param.angle   = [-angle_0/180*pi -angle_0/180*pi -angle_0/180*pi]; % radians
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_-90_-90_-90_trans_0_0_0.mat'};
    register(param)
    %close all;
    
    keyboard
    
    param.angle   = [angle_0/180*pi 0 angle_0/180*pi]; % radians
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_90_0_90_trans_0_0_0.mat'};
    register(param)
    close all;
    
    param.angle   = [angle_0/180*pi angle_0/180*pi angle_0/180*pi]; % radians
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_90_90_90_trans_0_0_0.mat'};
    register(param)
    close all;
    
    param.angle   = [0 angle_0/180*pi 0]; % radians
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_90_0_trans_0_0_0.mat'};
    register(param)
    close all;
    
    param.angle   = [0 0 angle_0/180*pi]; % radians
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_0_90_trans_0_0_0.mat'};
    register(param)
    close all;
    
    param.angle   = [angle_0/180*pi 0 0]; % radians
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_90_0_0_trans_0_0_0.mat'};
    register(param)
    close all;
    
    param.angle   = [0 0 0]; % radians
    
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_0_0_trans_0_0_50.mat'};
    register(param)
    close all;
    
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_0_0_trans_0_50_0.mat'};
    register(param)
    close all;
    
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_0_0_trans_50_0_0.mat'};
    register(param)
    close all;
    
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_0_0_trans_50_50_0.mat'};
    register(param)
    close all;
    
    param.inputFileName = {'bead_phantom_150_20180614_083753_rot_0_0_0_trans_50_50_50.mat'};
    register(param)
    close all;
end

%%
% for '/registration/repeat_90_90_0_D/', each repet reaches max MI. Does
% final rot and trans yield max MI for just_MI?

param.savePath = [param.inputFilePath1 '/registration/eval_repeat_90_90_0_D/'];

param.just_MI = true;
param.inputFileName = {'bead_phantom_150_20180614_083753_rot_90_90_0_trans_0_0_0.mat'};

param.angle = [1.5576 1.5757 -0.0166];
param.trans = [75.0171 75.2688 75.0222];
register(param)
close all;

param.angle = [1.5434 1.5662 -0.0266];
param.trans = [74.6410 75 75.1926];
register(param)
close all;

param.angle = [1.5389 1.5613 -0.0259];
param.trans = [74.5689 75 75];
register(param)
close all;

param.angle = [1.5576 1.5757 -0.0166];
param.trans = [75.0171 75.2688 75.0222];
register(param)
close all;

param.angle = [1.5434 1.5662 -0.0266];
param.trans = [74.6410 75 75.1926];
register(param)
close all;

param.angle = [1.5389 1.5613 -0.0259];
param.trans = [74.5689 75 75];
register(param)
close all;

param.angle = [1.5560 1.5699 -0.0140];
param.trans = [74.8604 75.0033 75];
register(param)
close all;

param.angle = [1.5421 1.5628 -0.0232];
param.trans = [74.6439 74.8713 75.0442];
register(param)
close all;

param.angle = [1.5798 1.5685 0.0096];
param.trans = [75.0358 75 75];
register(param)
close all;

param.angle = [1.5590 1.5653 -0.0091];
param.trans = [74.7200 75.0395 74.7720];
register(param)
close all;

param.angle = [1.5584 1.5758 -0.0145];
param.trans = [75.4586 74.8533 74.8612];
register(param)
close all;

param.angle = [1.5651 1.5714 -0.0011];
param.trans = [74.8041 74.9738 74.7073];
register(param)
close all;

param.angle = [1.5579 1.5634 -0.0095];
param.trans = [74.6423 74.7215 74.9853];
register(param)
close all;
