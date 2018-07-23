function [d,r] = perturb (param, gain)
d1 = random('unif',-1,1) * param.trans_amp;
d2 = random('unif',-1,1) * param.trans_amp;
d3 = random('unif',-1,1) * param.trans_amp;
r1 = random('unif',-1,1) * param.rot_amp(1) * gain;
r2 = random('unif',-1,1) * param.rot_amp(2) * gain;
r3 = random('unif',-1,1) * param.rot_amp(3) * gain;
d = [d1 d2 d3];
r = [r1 r2 r3];
%r = [r1 0 0];
end