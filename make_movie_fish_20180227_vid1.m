clear all;
close all;

%% load parameters
param = init_param();

%% setup parameters
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

if ~isempty(strfind(param.hostname, 'Justins-Mac'))
        param.ppath = '/Users/justin/Desktop/DLFM';
        addpath([param.ppath '/lotus-registration']);
        param.ipath = [param.ppath '/fish/20180227/vid1'];
        param.opath = param.ipath;
        param.spath = param.opath;
        param.inter = [param.ipath '/interpolate/'];
elseif ~isempty(strfind(param.hostname, 'willis'))
        param.ppath = '/home/jkinney/Desktop/DLFM';
        addpath([param.ppath '/lotus-registration']);
        param.ipath = [param.ppath '/fish/20180227/vid1'];
        param.opath = param.ipath;
        param.spath = param.opath;
        param.inter = [param.ipath '/interpolate/'];
else
        param.ppath = '/home/jkinney';
        addpath([param.ppath '/lotus-registration']);
 	param.ipath = '/om/project/boyden/DualLensLightField/01_23_18_NLS/video1';
        param.opath = '/om/scratch/Tue/jkinney/fish/20180227/vid1';
        param.spath = '/om/user/jkinney/DLFM/fish/20180227/vid1';
        param.inter = [param.spath '/interpolate/'];
end

param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];

param.savePath = [param.opath '/registration/'];

param.m = 1;
param.n = 500;

param.fpattern = 'Recon3D_solver_1_FrameNumber_%04d.mat';

param

prefix = 'fish_20180227_vid1';
timestamp = datestr(datetime('now'),'yyyymmdd_HHMMSS');
outputFileName = [param.spath '/' prefix '_' timestamp]

genVideo(param, outputFileName)
