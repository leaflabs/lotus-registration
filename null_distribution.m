function [cdf, centers, nullMIvec, param] = null_distribution (LFM1, LFM2, canonical, param)
nullf = [param.savePath param.inputFileName{1}(1:end-4) '_null.mat'];
if exist(nullf,'file') == 2
    load(nullf,'nullMIvec','bestMI','bestd','bestr');
    param.bestMI = bestMI;
    param.bestd = bestd;
    param.bestr = bestr;
else
    nullMIvec = [];
    param.bestMI = 0;
    param.bestd = [];
    param.bestr = [];
end
LL = length(nullMIvec);
fprintf('\nnull distribution has N = %d\n',LL);
% determine gain so offset limit = half or quarter of volume limits
% and so that rotation limit = pi
a = max ([ size(LFM1) size(LFM2)]);
offset_limit = 0.15 * a * param.voxel_y;
tmpt = param.trans_amp;
param.trans_amp = offset_limit;
%tmp = param.rot_amp;
gain = pi;
tmpr = param.rot_amp;
if param.confocal
    param.rot_amp = [1 1 1];
else
    param.rot_amp = [1 0 0];
end
%profile on;
fprintf('\nCount    Mutual_Information                Offset [um]                        Rotation [radians]\n');
% if param.parallel
%     parfor i=1:param.Nnull
%         MI = 0;
%         % perturb pos
%         % randomly pick a translation vector and rotation vector
%         % to be added to current location
%         [d,r] = perturb(param,gain);
%         % apply transformation
%         % rotate an amount r PLUS param.rot
%         % thus param.rot tracks the current rotation
%         %rotated = rotate (canonical,r);
%         % translate rotated by an amount param.trans+d
%         % thus param.trans tracks the current position
%         new = translate (rotate (canonical, r), param.centroid+d);
%         %new = translate (rotated, param.centroid+d);
%         % measure mutual_information
%         if strcmp(param.myfunc_MI,'multiply')
%             MI = mutual_information (LFM1, new, LFM2, param);
%             %     elseif strcmp(param.myfunc_MI,'multiply_sqrt')
%             %         MI = mutual_information_sqrt (LFM1, new, LFM2, param);
%         else
%             disp('WTF!');
%             keyboard;
%         end
%         nullMIvec = [nullMIvec MI];
%         str = sprintf('i = %4d, MI = %16.0f, d = [%7.3f0 %7.3f0 %7.3f0], r = [%7.3f0 %7.3f0 %7.3f0]',i,MI,d(1),d(2),d(3),r(1),r(2),r(3));
%         disp(str);
%         %if MI>bestMI
%         %    bestMI = MI;
%         %    bestd = d;
%         %    bestr = r;
%         %end
%         %         profile off
%         %         profile viewer
%     end
% else
    for i=1:param.Nnull
        MI = 0;
        % perturb pos
        % randomly pick a translation vector and rotation vector
        % to be added to current location
        [d,r] = perturb(param, gain, gain);
        % apply transformation
        % rotate an amount r PLUS param.rot
        % thus param.rot tracks the current rotation
        %rotated = rotate (canonical,r);
        % translate rotated by an amount param.trans+d
        % thus param.trans tracks the current position
        new = translate (rotate (canonical, r), param.centroid+d);
        %new = translate (rotated, param.centroid+d);
        % measure mutual_information
        if strcmp(param.myfunc_MI,'multiply')
            MI = mutual_information (LFM1, new, LFM2, param, 0);
            %     elseif strcmp(param.myfunc_MI,'multiply_sqrt')
            %         MI = mutual_information_sqrt (LFM1, new, LFM2, param);
        else
            disp('WTF!');
            keyboard;
        end
        nullMIvec = [nullMIvec MI];
        str = sprintf('i = %4d, MI = %16.0f, d = [%7.3f0 %7.3f0 %7.3f0], r = [%7.3f0 %7.3f0 %7.3f0]',i,MI,d(1),d(2),d(3),r(1),r(2),r(3));
        disp(str);
        if MI>param.bestMI
           param.bestMI = MI;
           param.bestd = d;
           param.bestr = r;
        end
        %         profile off
        %         profile viewer
    end
%end
param.rot_amp = tmpr;
param.trans_amp = tmpt;

% add to existing null distribution if any
fprintf('\nnull distribution has N = %d, was N = %d\n',length(nullMIvec),LL);
bestMI = param.bestMI;
bestd = param.bestd;
bestr = param.bestr;
save(nullf,'nullMIvec','bestMI','bestd','bestr');
fprintf('\nnull distribution saved as %s\n',nullf);

if ~param.parallel
    fprintf('\nBest of null:\n');
    fprintf('MI = %d\n',param.bestMI);
    fprintf('d = [%f %f %f]\n',param.bestd(1),param.bestd(2),param.bestd(3));
    fprintf('r = [%f %f %f]\n\n',param.bestr(1),param.bestr(2),param.bestr(3));
end

% CDF
centers = linspace(min(nullMIvec),max(nullMIvec),param.Nnull);
del = 0.5 * (centers(2)-centers(1));
edges  = [centers-del centers(end)+del];
%centers = 0.5 * (edges(1:end-1)+edges(2:end));
% this works. to check: sum(h_LFM1) == numel(LFM1)
h_MI = histcounts(nullMIvec,edges);
% normalize cdf to 1
h_MI = h_MI / sum(h_MI);

total = 0;
cdf = [];
for i=1:length(h_MI)
    total = total + h_MI(i);
    cdf = [cdf total];
end

% plot PDF and CDF
f = figure;
histogram(nullMIvec,50);
if strcmp(param.myfunc_MI,'multiply')
    xlabel('mutual information = sum(LFM1*LFM2)');
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    xlabel('mutual information = sum(sqrt(LFM1*LFM2))');
else
    disp('WTF!');
    keyboard;
end
ylabel('count');
title(['null distribution (bootstrapped from'...
    sprintf(' %d random registrations)',length(nullMIvec))]);

str=sprintf('%s%s_null_distribution.png',param.savePath,param.timestamp);
save_plot(f, str);

f = figure;
plot(centers,cdf);
if strcmp(param.myfunc_MI,'multiply')
    xlabel('mutual information = sum(LFM1*LFM2)');
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    xlabel('mutual information = sum(sqrt(LFM1*LFM2))');
else
    disp('WTF!');
    keyboard;
end
ylabel('normalized cumulative density function');
ylim([0 1]);
title(['null distribution (bootstrapped from'...
    sprintf(' %d random registrations)',length(nullMIvec))]);

str=sprintf('%s%s_null_cdf.png',param.savePath,param.timestamp);
save_plot(f, str);
end
