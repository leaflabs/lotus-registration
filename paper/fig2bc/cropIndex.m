function crop_ind = cropIndex(config, a, row, crop_row, z, crop_z)
dr = ceil(crop_row/2);
dz = ceil(crop_z/2);
rlim = [row-dr row+dr];
zlim = [z-dz z+dz];
% check bounds
rlim(1) = boundsCheck (rlim(1)<1, rlim(1), 1, 'Min row', a);
rlim(2) = boundsCheck (rlim(2)>a(1), rlim(2), a(1), 'Max row', a);
zlim(1) = boundsCheck (zlim(1)<1, zlim(1), 1, 'Min z', a);
zlim(2) = boundsCheck (zlim(2)>a(2), zlim(2), a(2), 'Max z', a);
if find([rlim(1) rlim(2) zlim(1) zlim(2)] < 0)
    crop_ind = -1;
    return
end
crop_ind = [rlim(1) rlim(2) zlim(1) zlim(2)];
end