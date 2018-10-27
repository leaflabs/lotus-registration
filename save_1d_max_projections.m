function param = save_1d_max_projections (LFM1, LFM2, new, param, str)
colors = {[0 0 1]
    [0.8 0.8 0.8]
    [0 0 1]
    [1 0 0]};

onetwo = squeeze(max(LFM1,[],3));
twothree = squeeze(max(LFM1,[],1));
LFM1_d1 = single(squeeze(max(onetwo,[],2)));
LFM1_d2 = single(squeeze(max(onetwo,[],1)));
LFM1_d3 = single(squeeze(max(twothree,[],1)));

onetwo = squeeze(max(LFM2,[],3));
twothree = squeeze(max(LFM2,[],1));
LFM2_d1 = single(squeeze(max(onetwo,[],2)));
LFM2_d2 = single(squeeze(max(onetwo,[],1)));
LFM2_d3 = single(squeeze(max(twothree,[],1)));

% projections for new
% convert new from microns to pixels
a = size(new);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_new = ceil(new.*scale);
ns = max(abc_new);

% create 2d projections of new in pixel space
onetwo = zeros(ns(1),ns(2));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,2) > 0
        onetwo(abc_new(i,1),abc_new(i,2)) = max([onetwo(abc_new(i,1),abc_new(i,2)) double(LFM2(param.index2(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(new)
    if abc_new(i,2) > 0 && abc_new(i,3) > 0
        twothree(abc_new(i,2),abc_new(i,3)) = max([twothree(abc_new(i,2),abc_new(i,3))  double(LFM2(param.index2(i)))]);
    end
end
% convert 2d projections to 1d
new_d1 = squeeze(max(onetwo,[],2));
new_d2 = squeeze(max(onetwo,[],1));
new_d3 = squeeze(max(twothree,[],1));

if 0 < 1
    T = '      (scaled by total intensity in sample)';
    scale_d1 = [ sum(LFM1_d1) sum(LFM2_d1) sum(new_d1)];
    scale_d2 = [ sum(LFM1_d2) sum(LFM2_d2) sum(new_d2)];
    scale_d3 = [ sum(LFM1_d3) sum(LFM2_d3) sum(new_d3)];
    yl = 'normalized intensity';
    %disp(T);
elseif 0 < 1
    T = '      (scaled by max intensity in each projected volume)';
    scale_d1 = single([ max(LFM1_d1) max(LFM2_d1) max(new_d1)]);
    scale_d2 = single([ max(LFM1_d2) max(LFM2_d2) max(new_d2)]);
    scale_d3 = single([ max(LFM1_d3) max(LFM2_d3) max(new_d3)]);
    yl = 'normalized intensity';
    %disp(T);
else
    T = '      (not scaled)';
    scale_d1 = [1 1 1 1 1];
    scale_d2 = [1 1 1 1 1];
    scale_d3 = [1 1 1 1 1];
    yl = 'intensity';
    %disp(T);
end

f = figure;
set(gcf,'Position',[79          18        1270         940]);
subplot(3,1,1);
semilogy(LFM1_d1/scale_d1(1),'.','Color',colors{1});
hold on;
semilogy(LFM2_d1/scale_d1(2),'.','Color',colors{2});
semilogy(new_d1/scale_d1(3),'-','Color',colors{4});
xlabel('first dimension [pixels]');
ylabel(yl);
hold off;
legend('LFM1','LFM2','LFM2 coarse reg');
title([str T],'Interpreter','none');
if param.xlim_1d(1) > 0
    xlim([0 param.xlim_1d(1)]);
else
    a = xlim;
    xlim([0 a(2)]);
    param.xlim_1d(1)=a(2);
end

subplot(3,1,2);
semilogy(LFM1_d2/scale_d2(1),'.','Color',colors{1});
hold on;
semilogy(LFM2_d2/scale_d2(2),'.','Color',colors{2});
semilogy(new_d2/scale_d2(3),'-','Color',colors{4});
xlabel('second dimension [pixels]');
ylabel(yl);
hold off;
legend('LFM1','LFM2','LFM2 coarse reg');
if param.xlim_1d(2) > 0
    xlim([0 param.xlim_1d(2)]);
else
    a = xlim;
    xlim([0 a(2)]);
    param.xlim_1d(2)=a(2);
end

subplot(3,1,3);
semilogy(LFM1_d3/scale_d3(1),'.','Color',colors{1});
hold on;
semilogy(LFM2_d3/scale_d3(2),'.','Color',colors{2});
semilogy(new_d3/scale_d3(3),'-','Color',colors{4});
xlabel('third dimension [pixels]');
ylabel(yl);
hold off;
legend('LFM1','LFM2','LFM2 coarse reg');
if param.xlim_1d(3) > 0
    xlim([0 param.xlim_1d(3)]);
else
    a = xlim;
    xlim([0 a(2)]);
    param.xlim_1d(3)=a(2);
end

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
save_plot(f, str);
end
