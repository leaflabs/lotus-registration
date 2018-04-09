clear all;
close all;

f = '/Users/justin/Desktop/DDLFM/fish/20180227/confocal_data/tiff_stack_video1_and_2.tif';
info = imfinfo(f);
num_images = numel(info);
width = info(1).Width;
height = info(1).Height;
data = zeros(height, width, num_images,'uint16');
for i=1:num_images
    fprintf('%d of %d\n',i,num_images);
    data(:,:,i) = uint16(imread(f,i));
end

% since DDLFM has voxel size
%
% param.voxel_x = 0.323; % um
% param.voxel_y = 0.323; % um
% param.voxel_z = 0.5;     % um
%
% and confocal has voxel size
%
% param.voxel_x_cf = 0.65; % um
% param.voxel_y_cf = 0.65; % um
% param.voxel_z_cf = 2.0;     % um
%
% interpolate confocal
% 2 in x
% 2 in y
% 4 in z
disp('Beginning interpolation 1');
tmp = interpolate (data, 2, 1);
disp('Beginning interpolation 2');
data = interpolate (tmp, 2, 2);
clear tmp;
disp('Beginning interpolation 3');
Xvolume = interpolate (data, 4, 3);
clear data;
disp('Saving data');
fout = [f(1:end-4) '.mat'];
save(fout, 'Xvolume', '-v7.3');
fprintf('# Output file = %s\n',fout);

function out = interpolate (vol, scale, dim)
% initialize container of new size
s = size(vol);
if dim == 1
    out = zeros(s(1)*scale,s(2),s(3),'uint16');
elseif dim == 2
    out = zeros(s(1),s(2)*scale,s(3),'uint16');
elseif dim == 3
    out = zeros(s(1),s(2),s(3)*scale,'uint16');
else
    disp('WTF?');
    keyboard
end
boundary = scale/2 + 1;
for i=1:scale*s(dim)
    if i < boundary
        if dim == 1
            out(i,:,:) = vol(1,:,:);
        elseif dim == 2
            out(:,i,:) = vol(:,1,:);
        elseif dim == 3
            out(:,:,i) = vol(:,:,1);
        else
            disp('WTF?');
            keyboard
        end
    elseif i > (scale*s(dim)-(boundary-1))
        if dim == 1
            out(i,:,:) = vol(end,:,:);
        elseif dim == 2
            out(:,i,:) = vol(:,end,:);
        elseif dim == 3
            out(:,:,i) = vol(:,:,end);
        else
            disp('WTF?');
            keyboard
        end
    else
        a = round((i-1)/scale);
        b = a+1;
        N = 2*scale;
        fb = 1+2*mod(i-(scale-1),scale);
        fa = N-fb;
        if dim == 1
            out(i,:,:) = fa/N*vol(a,:,:) + fb/N*vol(b,:,:);
        elseif dim == 2
            out(:,i,:) = fa/N*vol(:,a,:) + fb/N*vol(:,b,:);
        elseif dim == 3
            out(:,:,i) = fa/N*vol(:,:,a) + fb/N*vol(:,:,b);
        else
            disp('WTF?');
            keyboard
        end
    end
end
end



