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
        param.ipath = [param.ppath '/worm/20171215'];
        param.opath = param.ipath;
        param.inter = [param.ipath '/interpolate/'];
elseif ~isempty(strfind(param.hostname, 'willis'))
        param.ppath = '/home/jkinney/Desktop/DLFM';
        addpath([param.ppath '/lotus-registration']);
        param.ipath = [param.ppath '/worm/20171215'];
        param.opath = param.ipath;
        param.inter = [param.ipath '/interpolate/'];
else
        param.ppath = '/om/user/jkinney/DLFM';
        addpath([param.ppath '/lotus-registration']);
        param.ipath = '/om/project/boyden/DualLensLightField/12_15_17_40X_worm/video 1';
        param.opath = '/om/scratch/jkinney/worm/20171215';
        param.inter = [param.ipath '/interpolate/'];
end
param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

param.savePath = [param.opath '/registration/'];
param.voxel_x = 0.323/2; % um
param.voxel_y = 0.323/2; % um
param.voxel_z = 2.0;     % um

m = 1;
n = 1000;

param.rapid = true;

%% read log files
pathNpattern = sprintf('%s*FrameNumber_%04d_*.log',param.inter,m);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end
content = fileread([logs(1).folder '/' logs(1).name]);

eval(['transM = ' get_match('trans',content)]);
eval([  'rotM = ' get_match('rot',content)]);
eval([  'centroidM = ' get_match('centroid',content)]);
eval([  'clipM = ' get_match('clip',content)]);

pathNpattern = sprintf('%s*FrameNumber_%04d_*.log',param.inter,n);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end
content = fileread([logs(1).folder '/' logs(1).name]);

eval(['transN = ' get_match('trans',content)]);
eval([  'rotN = ' get_match('rot',content)]);
eval([  'centroidN = ' get_match('centroid',content)]);
eval([  'clipN = ' get_match('clip',content)]);

if ~isempty(find( (clipN-clipM) ~= 0))
    disp('WTF?');
    keyboard
end
param.clip = clipN; % equal to clipM
if ~isempty(find( (centroidN-centroidM) ~= 0))
    disp('WTF?');
    keyboard
end
param.centroid = centroidN ...
        + [param.clip(1)*param.voxel_y param.clip(3)*param.voxel_x 0];

for i=[param.m:param.n]
    param.rot = rotM + (rotN-rotM) * (i-m)/(n-m);
    tmp_trans = transM + (transN-transM) * (i-m)/(n-m);
    param.trans = tmp_trans + [ param.clip(1)*param.voxel_y  ...
        param.clip(3)*param.voxel_x   param.clip(5)*param.voxel_z ];
    param.inputFileName = {sprintf('Recon3D_solver_1_FrameNumber_%04d.mat',i)};
    
    param
    
    register(param)
end

function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{1}{1};
end


