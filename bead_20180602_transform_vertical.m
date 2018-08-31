clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
ppath = '/Users/justin/Desktop/DDLFM/bead/20180602';

savePath = [ppath '/vertical_registered'];
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

param.voxel_x = 0.5; % um
param.voxel_y = 0.5; % um
param.voxel_z = 4.0;     % um

p = [ppath '/vertical/Reconstructed/*.mat'];
d = dir(p);
fprintf('%d mat files found for\n%s\n\n',numel(d),p);

% from DDLFM/bead/20180602/summary.txt
%
% registration/Recon3D_after_1_Mono_N15__20180626_014848.log:                     rot: [1.5666 0.0027 -0.0207]
% registration/Recon3D_after_2_Mono_N15__20180702_212133.log:                     rot: [1.5789 0.0021 -0.0244]
% registration/Recon3D_after_3_Mono_N15__20180703_184504.log:                     rot: [1.5396 0.0123 -0.0132]
% registration/Recon3D_after_4_Mono_N15__20180704_041045.log:                     rot: [1.5882 0.0045 -0.0132]
% registration/Recon3D_after_5_Mono_N15__20180704_124421.log:                     rot: [1.5237 -0.0026 -0.0252]
% registration/Recon3D_after_6_Mono_N15__20180704_182133.log:                     rot: [1.4998 0.0072 -0.0170]
% registration/Recon3D_after_7_Mono_N15__20180705_015344.log:                     rot: [1.5356 -0.0051 -0.0363]
% registration/Recon3D_after_8_Mono_N15__20180705_124836.log:                     rot: [1.5189 -0.0021 -0.0195]
% registration/Recon3D_after_9_Mono_N15__20180708_171453.log:                     rot: [1.5745 0.0152 -0.0199]
% registration/Recon3D_after_10_Mono_N15__20180709_183438.log:                     rot: [1.5331 0.0137 -0.0323]
% registration/Recon3D_before_1_Mono_N15__20180621_112156.log:                     rot: [1.5534 0.0120 0.0086]
% registration/Recon3D_before_2_Mono_N15__20180629_190018.log:                     rot: [1.5625 0.0326 -0.0110]
% registration/Recon3D_before_3_Mono_N15__20180703_085911.log:                     rot: [1.6175 -0.0067 -0.0140]
% registration/Recon3D_before_4_Mono_N15__20180704_034338.log:                     rot: [1.6078 0.0206 9.8707e-05]
% registration/Recon3D_before_5_Mono_N15__20180704_065330.log:                     rot: [1.5223 0.0228 -0.0190]
% registration/Recon3D_before_6_Mono_N15__20180704_143222.log:                     rot: [1.5868 0.0223 -0.0123]
% registration/Recon3D_before_7_Mono_N15__20180704_184535.log:                     rot: [1.5489 0.0219 -0.0174]
% registration/Recon3D_before_8_Mono_N15__20180705_102751.log:                     rot: [1.6052 0.0069 -0.0043]
% registration/Recon3D_before_9_Mono_N15__20180705_145626.log:                     rot: [1.5376 0.0270 -0.0157]
% registration/Recon3D_before_10_Mono_N15__20180709_145804.log:                     rot: [1.5471 0.0128 -0.0095]
% 
% registration/Recon3D_after_1_Mono_N15__20180626_014848.log:                   trans: [259.9556 187.5921 182.7475]
% registration/Recon3D_after_2_Mono_N15__20180702_212133.log:                   trans: [259.5155 179.9928 188.7824]
% registration/Recon3D_after_3_Mono_N15__20180703_184504.log:                   trans: [258.9674 175.6954 184.4167]
% registration/Recon3D_after_4_Mono_N15__20180704_041045.log:                   trans: [260.6057 181.5537 185.4625]
% registration/Recon3D_after_5_Mono_N15__20180704_124421.log:                   trans: [258.6523 185.9795 199.2510]
% registration/Recon3D_after_6_Mono_N15__20180704_182133.log:                   trans: [259.4745 181.7909 191.1606]
% registration/Recon3D_after_7_Mono_N15__20180705_015344.log:                   trans: [257.8027 182.7796 194.6259]
% registration/Recon3D_after_8_Mono_N15__20180705_124836.log:                   trans: [258.4585 182.0175 191.7975]
% registration/Recon3D_after_9_Mono_N15__20180708_171453.log:                   trans: [259.7683 180.6478 186.1210]
% registration/Recon3D_after_10_Mono_N15__20180709_183438.log:                   trans: [260.5165 176.3722 193.7105]
% registration/Recon3D_before_1_Mono_N15__20180621_112156.log:                   trans: [258.2451 189.8673 195.7486]
% registration/Recon3D_before_2_Mono_N15__20180629_190018.log:                   trans: [261.0674 185.9299 190.3733]
% registration/Recon3D_before_3_Mono_N15__20180703_085911.log:                   trans: [257.3322 180.8416 196.0864]
% registration/Recon3D_before_4_Mono_N15__20180704_034338.log:                   trans: [258.9413 177.6604 192.9604]
% registration/Recon3D_before_5_Mono_N15__20180704_065330.log:                   trans: [260.7343 189.2274 190.4853]
% registration/Recon3D_before_6_Mono_N15__20180704_143222.log:                   trans: [260.9314 174.3025 188.3498]
% registration/Recon3D_before_7_Mono_N15__20180704_184535.log:                   trans: [260.6411 180.6617 192.0295]
% registration/Recon3D_before_8_Mono_N15__20180705_102751.log:                   trans: [258.9472 189.2522 191.5023]
% registration/Recon3D_before_9_Mono_N15__20180705_145626.log:                   trans: [260.4303 187.5672 190.3679]
% registration/Recon3D_before_10_Mono_N15__20180709_145804.log:                   trans: [259.6898 180.6363 192.5790]

d = {};
% d(1).name = 'Recon3D_after_1_Mono_N15.mat';
% d(1).rot = [1.5666 0.0027 -0.0207];
% d(1).trans = [259.9556 187.5921 182.7475];
d(1).name = 'Recon3D_after_1_Mono_N15.mat';
d(1).rot =  [1.5666 0.0027 -0.0207];
d(1).trans =  [259.9556 187.5921 182.7475];
d(2).name = 'Recon3D_after_2_Mono_N15.mat';
d(2).rot =  [1.5789 0.0021 -0.0244];
d(2).trans =  [259.5155 179.9928 188.7824];
d(3).name = 'Recon3D_after_3_Mono_N15.mat';
d(3).rot =  [1.5396 0.0123 -0.0132];
d(3).trans =  [258.9674 175.6954 184.4167];
d(4).name = 'Recon3D_after_4_Mono_N15.mat';
d(4).rot =  [1.5882 0.0045 -0.0132];
d(4).trans =  [260.6057 181.5537 185.4625];
d(5).name = 'Recon3D_after_5_Mono_N15.mat';
d(5).rot =  [1.5237 -0.0026 -0.0252];
d(5).trans =  [258.6523 185.9795 199.2510];
d(6).name = 'Recon3D_after_6_Mono_N15.mat';
d(6).rot =  [1.4998 0.0072 -0.0170];
d(6).trans =  [259.4745 181.7909 191.1606];
d(7).name = 'Recon3D_after_7_Mono_N15.mat';
d(7).rot =  [1.5356 -0.0051 -0.0363];
d(7).trans =  [257.8027 182.7796 194.6259];
d(8).name = 'Recon3D_after_8_Mono_N15.mat';
d(8).rot =  [1.5189 -0.0021 -0.0195];
d(8).trans =  [258.4585 182.0175 191.7975];
d(9).name = 'Recon3D_after_9_Mono_N15.mat';
d(9).rot =  [1.5745 0.0152 -0.0199];
d(9).trans =  [259.7683 180.6478 186.1210];
d(10).name = 'Recon3D_after_10_Mono_N15.mat';
d(10).rot =  [1.5331 0.0137 -0.0323];
d(10).trans =  [260.5165 176.3722 193.7105];
d(11).name = 'Recon3D_before_1_Mono_N15.mat';
d(11).rot =  [1.5534 0.0120 0.0086];
d(11).trans =  [258.2451 189.8673 195.7486];
d(12).name = 'Recon3D_before_2_Mono_N15.mat';
d(12).rot =  [1.5625 0.0326 -0.0110];
d(12).trans =  [261.0674 185.9299 190.3733];
d(13).name = 'Recon3D_before_3_Mono_N15.mat';
d(13).rot =  [1.6175 -0.0067 -0.0140];
d(13).trans =  [257.3322 180.8416 196.0864];
d(14).name = 'Recon3D_before_4_Mono_N15.mat';
d(14).rot =  [1.6078 0.0206 9.8707e-05];
d(14).trans =  [258.9413 177.6604 192.9604];
d(15).name = 'Recon3D_before_5_Mono_N15.mat';
d(15).rot =  [1.5223 0.0228 -0.0190];
d(15).trans =  [260.7343 189.2274 190.4853];
d(16).name = 'Recon3D_before_6_Mono_N15.mat';
d(16).rot =  [1.5868 0.0223 -0.0123];
d(16).trans =  [260.9314 174.3025 188.3498];
d(17).name = 'Recon3D_before_7_Mono_N15.mat';
d(17).rot =  [1.5489 0.0219 -0.0174];
d(17).trans =  [260.6411 180.6617 192.0295];
d(18).name = 'Recon3D_before_8_Mono_N15.mat';
d(18).rot =  [1.6052 0.0069 -0.0043];
d(18).trans =  [258.9472 189.2522 191.5023];
d(19).name = 'Recon3D_before_9_Mono_N15.mat';
d(19).rot =  [1.5376 0.0270 -0.0157];
d(19).trans =  [260.4303 187.5672 190.3679];
d(20).name = 'Recon3D_before_10_Mono_N15.mat';
d(20).rot =  [1.5471 0.0128 -0.0095];
d(20).trans =  [259.6898 180.6363 192.5790];

clip_microns = [param.clip(1) param.clip(3) param.clip(5)].*[param.voxel_y param.voxel_x param.voxel_z];

% skip i==4

for i=6:20
    f = [ppath '/vertical/Reconstructed/' d(i).name];
    fprintf('Loading %s...\n',f);
    load(f,'XguessSAVE1');
    

    param.centroid = [268.7500 201.2500 202] + clip_microns;
    if 0 > 1
        disp('CHECK THIS')
        s = size(XguessSAVE1);
        out = s/2.*[param.voxel_y param.voxel_x param.voxel_z];
        keyboard
    end
    param.trans = d(i).trans + clip_microns;
    param.rot = d(i).rot;
    prefix = d(i).name;
    mystr = sprintf('%s/%s_%s',savePath,prefix(1:end-4),timestamp);
    out = interpolate (XguessSAVE1,param);
    param.voxel_z = param.voxel_z / param.interp;
    clear XguessSAVE1;
    param.i = i;
    param.d = d;
    transform(out, param, mystr);
end

diary off;

%%
function transform (A, param, mystr)
a = size(A);
fprintf('size of input volume = [%d %d %d]\n',a(1),a(2),a(3));
SHIFT = 200;
EXPAND = 300;
%out = zeros(a);
out = zeros(a+EXPAND,class(A));
a = size(out);
fprintf('size of output volume = [%d %d %d]\n',a(1),a(2),a(3));


E = 100;
F = numel(A)/E;
for i=1:E
    fprintf('loop %d of %d\n',i,E);
    myrange = [1+(i-1)*F : i*F ]';
    %pos = init_pos([1:numel(A)]', A, param);
    pos = init_pos(myrange, A, param);
    %param.centroid = calc_centroid(A,param);
    tmp = translate (pos, -param.centroid);
    tmp = rotate (tmp, param.rot); % CHECK
    %tmp = translate (tmp, param.centroid);
    pos = translate (tmp, param.trans);
    B = render (pos, param, SHIFT);
    if any(max(B)>a)
        fprintf('Error! Rendered max (volume B) is outside of container out.\n');
        zz = param.d;
        fprintf('Rerun %s\n',zz(param.i).name);
        return;
    end
    if any(min(B)<1)
        fprintf('Error! Rendered min (volume B) is outside of container out.\n');
        zz = param.d;
        fprintf('Rerun %s\n',zz(param.i).name);
        return;
    end
    %out(B) = A(myrange);
    for j=1:length(myrange)
        out(B(j,1),B(j,2),B(j,3)) = A(myrange(j));
    end
    %keyboard
end




%out(1,1,1) = 0.5;
%ind = sub2ind(size(LFM1),abc_LFM1(:,1),abc_LFM1(:,2),abc_LFM1(:,3));
%i1 = LFM1(ind);
%out = uint16( single(i1) );
% save as MAT file
outFile = [mystr '.mat'];
% outFile = sprintf('%s_rot_%d_%d_%d_trans_%d_%d_%d.mat',...
%     mystr,...
%     param.rot(1)*180/pi, param.rot(2)*180/pi, param.rot(3)*180/pi,...
%     param.trans(1), param.trans(2), param.trans(3));
fprintf('Saving volume to %s.\n',outFile);
save(outFile,'out','-v7.3');
% save as TIF file
% outFile = sprintf('%s_rot_%d_%d_%d_trans_%d_%d_%d.tif',...
%     mystr,...
%     param.rot(1)*180/pi, param.rot(2)*180/pi, param.rot(3)*180/pi,...
%     param.trans(1), param.trans(2), param.trans(3));
outFile = [mystr '.tif'];
fprintf('Saving volume to %s.\n',outFile);
save_vol( out, outFile);
end

%%
function save_vol (A, outFile)
imwrite( squeeze(A(:,:,1)), outFile);
for k = 2:size(A,3)
    imwrite(squeeze(A(:,:,k)), outFile, 'WriteMode', 'append');
end
end

%%
function out = render (chunk, param, SHIFT)
%
% See function out = combine (LFM1, LFM2, chunk, linind, param)
% chunk = x,y,z centroids of voxels in LFM2
% for each centroid in chunk, calculate index of corresponding voxels in LFM1
L = length(chunk);
scale = [1/param.voxel_y*ones(L,1) ...
    1/param.voxel_x*ones(L,1) ...
    1/param.voxel_z*ones(L,1)];
out = ceil(chunk.*scale)+SHIFT;
if length(out) ~= length(unique(out,'rows'))
    disp('Uh oh');
    fprintf('length(out) = %d, length(unique(out,rows)) = %d\n',length(out),length(unique(out,'rows')));
end
b = max(out);
fprintf('max of indices (x, y, z) = [%d %d %d]\n',b(1),b(2),b(3));
b = min(out);
fprintf('min of indices (x, y, z) = [%d %d %d]\n',b(1),b(2),b(3));
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


function out = interpolate (LFM, param)
% assume param.interp = 1 or even
if ~(param.interp==1 || mod(param.interp,2)==0)
    fprintf('Error. param.interp is assumed to 1 or an even number (2,4,6,etc)\n.');
    fprintf('param.interp = %d\n',param.interp);
    exit
end
if param.interp>1
    % initialize container of new size
    s = size(LFM);
    out = zeros(s(1),s(2),param.interp*s(3),class(LFM));
    boundary = param.interp/2 + 2;
    a = 1;
    b = 0;
    j = 1;
    A = LFM(:,:,j);
    B = LFM(:,:,j+1);
    del = 1/param.interp;
    N = param.interp*s(3)+1;
    last_v = a*A+b*B;
    for i=1:N
        if i < boundary
            if i>1
                out(:,:,i-1) = LFM(:,:,1);
            end
        elseif i > (N-(boundary-1))
            out(:,:,i-1) = LFM(:,:,end);
        else
            v = a*A+b*B;
            out(:,:,i-1) = (last_v + v)/2;
            last_v = v;
            a = a-del;
            b = b+del;
            if a==0
                a=1;
                b=0;
                j=j+1;
                A = LFM(:,:,j);
                B = LFM(:,:,j+1);
            end
        end
    end
else
    out = LFM;
end
end

