function out = combine (LFM1, LFM2, chunk, linind, param)
% chunk = x,y,z centroids of voxels in LFM2
% for each centroid in chunk, calculate index of corresponding voxels in LFM1
L = length(chunk);
scale = [1/param.voxel_y*ones(L,1) ...
    1/param.voxel_x*ones(L,1) ...
    1/param.voxel_z*ones(L,1)];
abc_LFM1 = ceil(chunk.*scale);
% purge any indices that lie outLFM of LFM1
ind = find(abc_LFM1<1);
[r1,~] = ind2sub(size(abc_LFM1),ind);
s = size(LFM1);
ind = find(abc_LFM1(:,1)>s(1));
[r2,~] = ind2sub(size(abc_LFM1),ind);
ind = find(abc_LFM1(:,2)>s(2));
[r3,~] = ind2sub(size(abc_LFM1),ind);
ind = find(abc_LFM1(:,3)>s(3));
[r4,~] = ind2sub(size(abc_LFM1),ind);
rem = unique([r1;r2;r3;r4]);
N = length(abc_LFM1);
f_not_used = length(rem)/N;
fprintf('[Voxels] Number %d, fraction used %1.3f, with negative index %1.3f, outLFM of LFM1 %1.3f\n',...
    N,1-f_not_used,length(r1)/N,(length(r2)+length(r3)+length(r4))/N);
% chunk is a systematic sweep in z dimensino of LFM2
% in contrast, abc_LFM1 are the voxels in LFM1 coordinate system of a
% rotated and translated LFM2, so not well behaved
% thus, compute intensity in LFM2 coordinate system
out = ones(length(chunk),1,'uint16');
%out = zeros(length(chunk),1,'uint16');
%
% how to handle rem entries
% in abc, replace with 1,1,1
abc_LFM1(rem,:) = [ones(length(rem),1) ones(length(rem),1) ones(length(rem),1)];
%
% in intensity calculation, set all rem entries to zero
ind = sub2ind(size(LFM1),abc_LFM1(:,1),abc_LFM1(:,2),abc_LFM1(:,3));
i1 = LFM1(ind);
i1(rem) = 0.0;
i2 = LFM2(linind(1):linind(2))';

if strcmp(param.myfunc_combine,'multiply')
    % multiply
    out = uint16( single(i1).*single(i2) );
elseif strcmp(param.myfunc_combine,'multiply_sqrt')
    % multiply
    out = uint16( sqrt(single(i1).*single(i2)) );
elseif strcmp(param.myfunc_combine,'i1')
    out = i1;
elseif strcmp(param.myfunc_combine,'i2')
    out = i2;
elseif strcmp(param.myfunc_combine,'sum')
    out = i1 + i2;
elseif strcmp(param.myfunc_combine,'min')
    % OR min OR probability
    out = uint16( min(i1,i2));
else
    disp('WTF?!');
    keyboard
end
s = size(LFM2);
out = reshape(out,s(1),s(2));
%max(max(max(out)))
end