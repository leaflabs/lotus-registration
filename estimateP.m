function p = estimateP (param, canonical, LFM1, LFM2, new, T, gain)
init_MI = mutual_information (LFM1, new, LFM2, param, 0);
Pvec = [];
%N = 100;
N = 10 / param.Pmelt;
j = N;
hot = 1;
while j>0
    % perturb pos
    % randomly pick a translation vector and rotation vector
    % to be added to current location
    [d,r] = perturb(param, gain, gain);
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
    %str0 = sprintf('d = [%7.3g0 %7.3g0 %7.3g0], r = [%7.3g0 %7.3g0 %7.3g0]',d(1),d(2),d(3),r(1),r(2),r(3));
    % apply transformation
    % rotate an amount r PLUS param.rot
    % thus param.rot tracks the current rotation
    rotated = rotate (canonical,param.rot+r);
    % translate rotated by an amount param.trans+d
    % thus param.trans tracks the current position
    new = translate (rotated, param.trans+d);
    %str1 = print_param(param);
    % measure mutual_information
    if strcmp(param.myfunc_MI,'multiply')
        MI = mutual_information (LFM1, new, LFM2, param, 0);
    else
        disp('WTF!');
        keyboard;
    end
    delmi = MI - init_MI;
    %str2 = sprintf('test MI = %7.3g, delmi = %7.3g',MI,delmi);
    %str3 = 'MI increased';
    if delmi < 0.0
        j=j-1;
        p = exp(delmi/T);
        rnd = rand(1);

        if rnd < p
            % accept decrease in MI mutual information
            %str3 = sprintf('Accept %.3g < %.3g',rnd,p);
            Pvec = [Pvec 1];
        else
            % reject move
            Pvec = [Pvec 0];
            %str3 = sprintf('Reject %.3g > %.3g',rnd,p);
        end
    end
    %str4 = sprintf('%d %d, T = %1.0e',i,j,T);
    %fprintf('%7s%75s  %84s  %40s  %22s\n',str4,str0,str1,str2,str3);
end
p = sum(Pvec)/N;
end
