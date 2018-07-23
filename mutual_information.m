function mi = mutual_information (LFM1, new, LFM2, param)
scale = [1/param.voxel_y 0 0 ; 0 1/param.voxel_x 0 ; 0 0 1/param.voxel_z];
new_pixels = ceil(new*scale);
s = size(LFM1);
% find rows of new_pixels that overlap LFM1 and save as index
index = find( new_pixels(:,1)<=s(1) & new_pixels(:,2)<=s(2) & new_pixels(:,3)<=s(3) ...
    & new_pixels(:,1)>0 & new_pixels(:,2)>0 & new_pixels(:,3)>0 );
% in LFM1, lookup intensity value at 
i1 = LFM1(sub2ind(s,new_pixels(index,1),new_pixels(index,2),new_pixels(index,3)));
i2 = LFM2(param.index2(index));
mi = sum(double(i1).*double(i2));
end