function [out] = rotate (pos, angle)
rot = single( rotation_matrix (angle) );
%out = pos*rot;
out = rot*pos';
out = out';
end