function T = find_melting_T (LFM1, new, LFM2, canonical, param, gain)
if param.T_fast
    T = param.T0;
    return;
end
%param.rot_amp = [param.rot_amp(1) 0 0];
% calculate MI
% if strcmp(param.myfunc_MI,'multiply')
%     MI = mutual_information (LFM1, new, LFM2, param);
% else
%     disp('WTF!');
%     keyboard;
% end
%init_MI = MI;
% set initial T
% Start with T sufficiently high to "melt" the system
% later to guarantee melting, if needed T will be increased until
% P (accepting a decrease in MI) = 60%
T = param.T0;
% while system not frozen and more temperature changes are allowed
% profile on;
%%%% TODO -Add condition 'if no change in voxel overlap between LFM1 + DDLFM
%lP = 0;
i = 0;
while 1 > 0
    i=i+1;
    fprintf('i = %d\n',i);
    p = estimateP (param, canonical, LFM1, LFM2, new, T, gain);
    fprintf('\nfraction accepted = %0.5f, T = %1.0e\n',p,T);
    if p < param.Pmelt
        T = T * 10;
        %lP = p;
        fprintf('new T = %1.0e\n\n',T);
    else
        break;
    end
end
hT = T;
lT = T/10;
mT = (hT+lT)/2;
i = 0;
mP = estimateP (param, canonical, LFM1, LFM2, new, mT, gain);
Pdiff = abs(mP-param.Pmelt);
fprintf('[hT = %1.5e, mT = %1.5e, lT = %1.5e, mP = %1.5f, Pmelt = %1.5f,Pdiff = %1.5f, Pepsilon = %1.5f\n',hT,mT,lT,mP,param.Pmelt,Pdiff,param.Pepsilon);
while Pdiff > param.Pepsilon
    i=i+1;
    if mP < param.Pmelt
        % increase temp
        lT = mT;
        if abs(lT-hT)<1e8
            hT = lT*10;
        end
    else
        % lower temp
        hT = mT;
        if abs(lT-hT)<1e8
            lT = hT/10;
        end
    end
    mT = (hT+lT)/2;
    mP = estimateP (param, canonical, LFM1, LFM2, new, mT, gain);
    Pdiff = abs(mP-param.Pmelt);
    fprintf('[hT = %1.5e, mT = %1.5e, lT = %1.5e, mP = %1.5f, Pmelt = %1.5f, Pdiff = %1.5f, Pepsilon = %1.5f\n',hT,mT,lT,mP,param.Pmelt,Pdiff,param.Pepsilon);
end
T = mT;
end
