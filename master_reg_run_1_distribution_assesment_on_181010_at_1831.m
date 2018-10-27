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
	
% elseif ~isempty(strfind(param.hostname, 'willis'))
% 	param.ppath = '/home/jkinney/Desktop/DDLFM';
% 	addpath([param.ppath '/lotus-registration']);
% 	param.ipath = [param.ppath '/bead/20180602'];
% 	param.opath = param.ipath;
%     param.parallel = false;
% else
% 	param.ppath = '/home/jkinney';
% 	addpath([param.ppath '/lotus-registration']);
% 	param.ipath = '/om/user/jkinney/DLFM/bead/20180602';
% 	param.opath = param.ipath;
%         %param.opath = '/om/scratch/Tue/jkinney/bead/20180602';
end

param.parallel = false;

%% 
param.voxel_x = 0.5; % um
param.voxel_y = 0.5; % um
param.voxel_z = 4;     % um

param.angle   = [-pi/2 0 0]; % radians

if 0 > 1
    %% Corban's settings
    % SHOULD BE param.interp = 1;
    param.Pmelt = 0.2;
    param.trans_amp = 3.5;
    param.rot_amp = [1 0.1 0.1];
    param.Tmin = 1e7;
    param.Nnull = 2000;
    param.T_fast = false;
    param.frozen2 = 100;
elseif 0 > 1
    %% all registrations Oct 23-24
    param.Pmelt = 0.3;
    param.Pdiff = 0.03;
    param.trans_amp = 3.5;
    param.rot_amp = [2 2 2];
    param.Tmin = 1e1;
    param.Nnull = 20;
    param.T_fast = true;
    param.frozen2 = 100;
    param.clip = zeros(1,6);
    param.dynamic_range_thresh = 2;
    param.T0 = 1e7;
    param.last_pass_N = 25;
    param.frozen3 = 6000;
    % last_pass_N should be less than frozen2, e.g. frozen2/2, I like 25
else
    %param.plot = false;
    %param.Pmelt = 0.3;
    %param.Pepsilon = 0.03;
    %param.trans_amp = 3.5;
    param.rot_amp = [2 2 2];
    %param.Tmin = 1e1;
    %param.Nnull = 20;
    %param.T_fast = false;
   % param.frozen2 = 100;
    param.clip = zeros(1,6);
    %param.dynamic_range_thresh = 2;
    %param.T0 = 1e7;
    %param.last_pass_N = 25;
    param.frozen3 = 6000;
    param.offset = [-13 -13 8];
    param.angle   = [-1.0*pi/2 0.03 -0.07];
    % last_pass_N should be less than frozen2, e.g. frozen2/2, I like 25
end

%for i=[12 15 24]
%for i=[2 5 9]

% phantom 2, data contains only three beads? are two beads missing?
% phantom 5, data contains only four beads? is one bead missing?
% phantom 6, try annealing for 10000 iterations 
% phantom 9, coarse trans in dim3 is way off
% phantom 13, coarse trans in dim2 is way off
% phantom 24, data contains only four beads? is one bead missing? 

%for i=[2 5 6 9 13 24]
for i=[7]
    partial_path = sprintf('/run_1_(distribution_assesment)_on_181010_at_1831/phantom_vol_%02d_process',i);
    param = prep_paths (param, partial_path);
    param
    register(param)
    close all;
end

function param = prep_paths (param, partial_path)
param.ipath = [param.ppath partial_path];
param.opath = param.ipath;
param.inputFilePath1 = [param.ipath '/view_1_for_registration/'];
param.inputFilePath2 = [param.ipath '/view_2_for_registration/'];
param.savePath = [param.opath '/registration/'];
param.inputFileName = {['recon_3d_im.mat']};
end
