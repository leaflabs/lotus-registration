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
elseif 0 > 1
    %% all registrations Oct 24-25
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
    
    for i=1:25
        partial_path = sprintf('/run_1_(distribution_assesment)_on_181010_at_1831/phantom_vol_%02d_process',i);
        param = prep_paths (param, partial_path);
        param
        register(param)
        close all;
    end
    
else
    %% all registrations Oct 27-29
    param.rot_amp = [2 2 2];
    param.clip = zeros(1,6);
    param.frozen3 = 6000;
    param.Nnull = 20;
    %param.offset = [-13 -13 8];
    %param.angle   = [-1.0*pi/2 0.03 -0.07];
end

% data sets to repeat
ivec = 1:13;

top_path = [param.ppath '/run_1_(distribution_assesment)_on_181010_at_1831'];
d = dir([top_path '/**/*.log']);

% for up to 10 seeds of random number generator
for j=1:3
    [median_offset , median_rot] = get_median( d, j:3:numel(d) );
    param.offset = median_offset;
    param.angle  = median_rot;
    rng(j);
    % for chosen data sets
    for i=ivec
        fprintf('RNG seed = %d\n',j);
        partial_path = sprintf('/run_1_(distribution_assesment)_on_181010_at_1831/phantom_vol_%02d_process',i);
        param = prep_paths (param, partial_path);
        param
        register(param)
        close all;
    end
end


function [median_offset , median_rot] = get_median ( d, vec )
offset = [];
rot = [];
for i=vec
    fpath = [d(i).folder '/' d(i).name];
    fprintf('Reading log file: %s\n',fpath);
    out = read_log (fpath);
    offset = [offset ; out.offsetF];
    rot = [rot ; out.rot];
end
median_offset = median(offset);
median_rot = median(rot);
end

%% read log file
function out = read_log (f)
out = [];

content = fileread(f);

eval(['out.trans = ' get_match('trans',content) ';']);
%eval(['out.angle = ' get_match('angle',content) ';']);
eval(['out.centroid = ' get_match('centroid',content) ';']);
%eval(['out.offsetI = ' get_match('offset',content) ';']);
out.offsetF = out.trans - out.centroid;
eval(['out.rot = ' get_match('rot',content) ';']);
end

function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{end}{1};
end

function param = prep_paths (param, partial_path)
param.ipath = [param.ppath partial_path];
param.opath = param.ipath;
param.inputFilePath1 = [param.ipath '/view_1_for_registration/'];
param.inputFilePath2 = [param.ipath '/view_2_for_registration/'];
param.savePath = [param.opath '/registration/'];
param.inputFileName = {['recon_3d_im.mat']};
end
