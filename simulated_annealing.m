function [new, param] = simulated_annealing (LFM1, new, LFM2, canonical, param)
param.transvec = [param.trans];
param.offsetvec = [param.offset];
param.rotvec = [param.rot];
param.Pvec = [];
%param.trans_amp = 3.5; %um, half diameter of neuron
% calc rotational gain
a = size(LFM2);
half_span = 0.5*[a(1)*param.voxel_y a(2)*param.voxel_x];
radius = sqrt( sum( half_span .* half_span ) );
gain = param.trans_amp / radius;
% calculate MI
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, new, LFM2, param);
else
    disp('WTF!');
    keyboard;
end
last_MI = MI;
param.MIvec = MI;
% CDF
param.cdfvec = [];
w = find(param.centers>last_MI,1); % keep only first instance
if ~isempty(w)
    param.cdfvec = [param.cdfvec param.cdf(w)];
else
    val = double(last_MI)/double(param.bestMI);
    param.cdfvec = [param.cdfvec val] ;
end
% set initial T
% Start with T sufficiently high to "melt" the system
% later to guarantee melting, if needed T will be increased until
% P (accepting a decrease in MI) = 60%
T = find_melting_T(LFM1, new, LFM2, canonical, param, gain);
param.Tvec = T;
%param = set_prate (MI,param);
%p = param.init_p;
% set max number of temperature changes and mean changes
Tchanges = param.TC0;
if strcmp(param.anneal,'exp')
    param.Trate = 10^(log10(1e-3)/Tchanges);
elseif  strcmp(param.anneal,'linear')
    param.Trate = T / Tchanges;
else
    disp('WTF!');
    keyboard;
end
fprintf('\n\nTrate = %f\n\n',param.Trate);

param.Dvec = [];
% while system not frozen and more temperature changes are allowed
% profile on;
disp('Count      Perturbation     Transformation      Overlap     Probability,Decision   ');
%%%% TODO -Add condition 'if no change in voxel overlap between LFM1 + DDLFM
i = 1;
pass = 0;
hot = 1;
phase = 2;
while 1 > 0
    if phase==2 && pass>param.frozen2
        fprintf('(\n\n%d passes in a row. The model is frozen.\n\n',param.frozen2);
        break;
    end
    if phase==1 && pass>param.frozen1
        fprintf('(\n\n%d passes in a row. Switching to fine search.\n\n',param.frozen1);
        phase = 2;
        pass = 0;
    end
    i = i+1;
    if strcmp(param.anneal,'exp')
        T = T * param.Trate;
    elseif  strcmp(param.anneal,'linear')
        T = T - param.Trate;
    else
        disp('WTF!');
        keyboard;
    end
%     if i > Tchanges
%         param.Tvec = [param.Tvec T];
%         break;
%     end
    if T < param.Tmin
        param.Tvec = [param.Tvec T];
        break;
    end
    % perturb pos
    % randomly pick a translation vector and rotation vector
    % to be added to current location
    [d,r] = perturb(param,gain);
    if phase == 2
        if hot<4
            r = zeros(1,3);
            t = r;
            t(hot) = d(hot);
            d = t;
        else
            d = zeros(1,3);
            t = d;
            t(hot-3) = r(hot-3);
            r = t;
        end
        hot = hot+1;
        if hot>6
            hot = 1;
        end
    end
    str0 = sprintf('d = [%7.3f %7.3f %7.3f], r = [%7.5f %7.5f %7.5f]',d(1),d(2),d(3),r(1),r(2),r(3));
    % apply transformation
    % rotate an amount r PLUS param.rot
    % thus param.rot tracks the current rotation
    rotated = rotate (canonical,param.rot+r);
    param.rot = param.rot + r;
    % translate rotated by an amount param.trans+d
    % thus param.trans tracks the current position
    new = translate (rotated, param.trans+d);
    param.trans = param.trans + d;
    str1 = print_param(param);
    % measure mutual_information
    if strcmp(param.myfunc_MI,'multiply')
        MI = mutual_information (LFM1, new, LFM2, param);
    else
        disp('WTF!');
        keyboard;
    end
    delmi = MI - param.MIvec(end);
    str2 = sprintf('test MI = %7.3g, delmi = %7.3g',MI,delmi);
    if delmi > 0.0
        % keep
        str3 = 'Accept - MI increased';
        param.MIvec = [param.MIvec MI];
        last_MI = MI;
        pass = 0;
    else
        %         # else accept D with probability P = exp(-E/kBT)
        %         # using (psuedo-)random number uniformly distributed in the interval (0,1)
        p = exp(delmi/T);
        param.Pvec = [param.Pvec p];
        % Specifically, if random number is less the P, then accept
        rnd = rand(1);
        if rnd < p
            % accept decrease in MI mutual information
            str3 = sprintf('Accept %.3g < %.3g',rnd,p);
            param.MIvec = [param.MIvec MI];
            last_MI = MI;
            if delmi==0
                pass = pass + 1;
            else
                pass = 0;
            end
        else
            % reject move
            str3 = sprintf('Reject %.3g > %.3g',rnd,p);
            param.rot = param.rot - r;
            param.trans = param.trans - d;
            tmp = param.MIvec(end);
            param.MIvec = [param.MIvec tmp];
            pass = pass + 1;
        end
    end
    str4 = sprintf('%d, pass = %d, T = %g',i,pass, T);
    str6 = sprintf('Final MI = %d',last_MI);
    str5 = print_param(param);
    val7 = param.trans - param.centroid;
    str7 = sprintf('final offset = [%6.6g0 %6.6g0 %6.6g0]',val7(1),val7(2),val7(3));
    w = find(param.centers>last_MI,1); % keep only first instance
    if last_MI >= param.bestMI
        val = double(last_MI)/double(param.bestMI);
        str8 = sprintf('MI frac = %.5g, siman exceeds null',val);
        param.cdfvec = [param.cdfvec val] ;
    elseif ~isempty(w)
        str8 = sprintf('MI frac = %.5g, null exceeds siman',param.cdf(w));
        param.cdfvec = [param.cdfvec param.cdf(w)];
    else
        disp('WTF?');
        keyboard
    end
    fprintf('%7s%75s  %84s  %40s  %22s  %84s  %20s %40s %20s\n',str4,str0,str1,str2,str3,str5,str6,str7,str8);
    param.Tvec = [param.Tvec T];
    param.transvec = [param.transvec; param.trans];
    param.offsetvec = [param.offsetvec; val7];
    param.rotvec = [param.rotvec; param.rot];
    param.Dvec = [param.Dvec; d];
end
rotated = rotate (canonical,param.rot);
new = translate (rotated, param.trans);
end