function pos = init_pos (linear_i,LFM, param)
[a,b,c] = ind2sub(size(LFM),linear_i);
y = single( (a-0.5) * param.voxel_y );
x = single( (b-0.5) * param.voxel_x );
z = single( (c-0.5) * param.voxel_z );
pos = [y x z];
end
