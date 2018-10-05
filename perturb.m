function [d,r] = perturb (param, gain, gain0)
% but gain is not used for d
% so introduce gain ratio to turn down d
a = gain / gain0;
d1 = random('unif',-1,1) * param.trans_amp * a;
d2 = random('unif',-1,1) * param.trans_amp * a;
d3 = random('unif',-1,1) * param.trans_amp * a;
% how gain was originally used
r1 = random('unif',-1,1) * param.rot_amp(1) * gain;
r2 = random('unif',-1,1) * param.rot_amp(2) * gain;
r3 = random('unif',-1,1) * param.rot_amp(3) * gain;
d = [d1 d2 d3];
r = [r1 r2 r3];
%r = [r1 0 0];
end
