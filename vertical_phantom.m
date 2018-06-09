clear all;
close all;

%% load parameters
param = init_param();

% without loss of generality
param.voxel_x = 1;
param.voxel_y = 1;
param.voxel_z = 1;

[~,hostname] = system('hostname')

if ~isempty(strfind(hostname, 'Justins-Mac'))
	ppath = '/Users/justin/Desktop/DDLFM/bead/20180328';
elseif ~isempty(strfind(hostname, 'willis'))
	ppath = '/home/jkinney/Desktop/DDLFM/bead/20180328';
else
	disp('WTF?');
    keyboard
end

savePath = [ppath '/vertical_phantom'];
if ~exist(savePath,'dir')
    status = mkdir(savePath);
    if status == 1
        disp(['Created folder: ' savePath]);
    else
        disp(['Error attempting to create folder:' savePath]);
        status
        exit;
    end
end

timestamp = [datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%s/output_%s.log',savePath, timestamp);
diary(fname)
tic

% p = [ppath '/*.mat'];
% d = dir(p);
% fprintf('%d mat files found for\n%s\n\n',numel(d),p);

%f = [d(1).folder '/' d(1).name];

postfix = 'Recon3D_after_9_Mono_N15.mat';

f = [ppath '/vertical/Reconstructed/' postfix];
fprintf('Loading %s...\n',f);
load(f,'XguessSAVE1');
A = XguessSAVE1;
clear XguessSAVE1;

mystr = sprintf('%s/%s_%s',savePath,postfix(1:end-4),timestamp);

% param.rot = [0 0 0];
% param.trans = [0 0 0];
% 
% transform(A, param, mystr);

param.rot = [pi/2 0 0];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [0 pi/2 0];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [0 0 pi/2];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [pi/2 pi/2 0];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [pi/2 0 pi/2];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [pi/2 pi/2 pi/2];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [-pi/2 -pi/2 -pi/2];
param.trans = [0 0 0];

transform(A, param, mystr);

param.rot = [0 0 0];
param.trans = [50 0 0];

transform(A, param, mystr);

param.rot = [0 0 0];
param.trans = [0 50 0];

transform(A, param, mystr);

param.rot = [0 0 0];
param.trans = [0 0 50];

transform(A, param, mystr);

param.rot = [0 0 0];
param.trans = [50 50 0];

transform(A, param, mystr);

param.rot = [0 0 0];
param.trans = [50 50 50];

transform(A, param, mystr);

%param.rot = [-pi/2 -pi/2 -pi/2 ];
%param.trans = [-5 -5 -5];

%transform(A, param, mystr);

diary off;

%%
function transform (A, param, mystr)
pos = init_pos([1:numel(A)]', A, param);
param.centroid = calc_centroid(A,param);
fprintf('Translating...\n');
tmp = translate (pos, -param.centroid);
fprintf('Rotating...\n');
tmp = rotate (tmp, param.rot); % CHECK
fprintf('Translating...\n');
tmp = translate (tmp, param.centroid+param.trans);
%fprintf('Translating...\n');
%pos = translate (tmp, param.trans);
A = render (A, pos, [1:numel(A)]', param);
% save as MAT file
outFile = sprintf('%s_rot_%d_%d_%d_trans_%d_%d_%d.mat',...
    mystr,...
    param.rot(1)*180/pi, param.rot(2)*180/pi, param.rot(3)*180/pi,...
    param.trans(1), param.trans(2), param.trans(3));
fprintf('Saving volume to %s.\n',outFile);
save(outFile,'A','-v7.3');
% save as TIF file
outFile = sprintf('%s_rot_%d_%d_%d_trans_%d_%d_%d.tif',...
    mystr,...
    param.rot(1)*180/pi, param.rot(2)*180/pi, param.rot(3)*180/pi,...
    param.trans(1), param.trans(2), param.trans(3));
fprintf('Saving volume to %s.\n',outFile);
save_vol( A, outFile);
end

%%
function save_vol (A, outFile)
imwrite( squeeze(A(:,:,1)), outFile);
for k = 2:size(A,3)
    imwrite(squeeze(A(:,:,k)), outFile, 'WriteMode', 'append');
end
end

%%
function out = render (LFM1, chunk, linind, param)
%
% See function out = combine (LFM1, LFM2, chunk, linind, param)
% chunk = x,y,z centroids of voxels in LFM2
% for each centroid in chunk, calculate index of corresponding voxels in LFM1
L = length(chunk);
scale = [1/param.voxel_y*ones(L,1) ...
    1/param.voxel_x*ones(L,1) ...
    1/param.voxel_z*ones(L,1)];
abc_LFM1 = ceil(chunk.*scale);
if length(abc_LFM1) ~= length(unique(abc_LFM1,'rows'))
    disp('Uh oh');
    keyboard
end

a = size(LFM1);
fprintf('size of input volume = [%d %d %d]\n',a(1),a(2),a(3));
b = max(abc_LFM1);
fprintf('max of indices (x, y, z) = [%d %d %d]\n',b(1),b(2),b(3));

%out = zeros(a);
out = zeros(b);

for i=1:L
    out(abc_LFM1(i,1),abc_LFM1(i,2),abc_LFM1(i,3)) = LFM1(i);
end
out(1,1,1) = 0.5;
%ind = sub2ind(size(LFM1),abc_LFM1(:,1),abc_LFM1(:,2),abc_LFM1(:,3));
%i1 = LFM1(ind);
%out = uint16( single(i1) );
end

%%
function out = calc_centroid (LFM, param)
s = size(LFM);
out = s/2.*[param.voxel_y param.voxel_x param.voxel_z];
end

%%
function [out] = rotate (pos, angle)
rot = single( rotation_matrix (angle) );
%out = pos*rot;
out = rot*pos';
out = out';
end

%%
function out = rotation_matrix (angle)
a = [1 0 0;...
    0 cos(angle(1)) -sin(angle(1));...
    0 sin(angle(1)) cos(angle(1))];
b = [cos(angle(2)) 0 sin(angle(2));...
    0 1 0;...
    -sin(angle(2)) 0 cos(angle(2))];
c = [cos(angle(3)) -sin(angle(3)) 0;...
    sin(angle(3)) cos(angle(3)) 0;...
    0 0 1];
%out = a*b*c;
out = c*b*a;
end

%%
function [out] = translate (pos, delta)
D = ones(size(pos),'single').*delta;
out = pos + D;
end

%%
function pos = init_pos (linear_i,LFM, param)
[a,b,c] = ind2sub(size(LFM),linear_i);
y = single( (a-0.5) * param.voxel_y );
x = single( (b-0.5) * param.voxel_x );
z = single( (c-0.5) * param.voxel_z );
pos = [y x z];
end