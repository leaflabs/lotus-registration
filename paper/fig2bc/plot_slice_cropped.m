function [dim,dim2,L2c] = plot_slice_cropped ( config, handle, ...
    csec_cropped, a, crop_index, xlbl, ylbl )
%%% render slice as image
subplot(handle);
bits = ceil(log2(single(max(max(csec_cropped)))));
cmap = gray(2^bits);
colormap(cmap);
im = imagesc(csec_cropped,[0 2^bits]);
daspect([1,1,1]);
%%% scale image according to pixel sizes
b = size(csec_cropped);
% vertical_microns = single(b(1)) * config.pixel;
% depth_microns = (single(b(2))-1.0) * single(config.zspacing);
% scale = vertical_microns / depth_microns;
% disp(['image size = ' num2str(vertical_microns) ' um tall x ' num2str(depth_microns) ' um wide' ]);
% pos = get(h,'Position'); % pos = [left bottom width height]
% newpos = pos;
% %newpos(1) = 10;
% newpos(3) = pos(3) * scale;
% set(h,'Position', newpos);
%%% relabel axes to microns
% to do so,
% specify XTick = tick locations in csec_cropped space
% specify XTickLabel = labels for each tick
%
% relative to cropped image (i.e. origin in upper left corner)
dim2 = [1];
L2 = [1];
dim1 = [1];
L1 = [1];
for i=1:config.div-1
    dim2 = [dim2 round(b(2)/config.div*i)];
    L2 = [L2 dim2(end)*config.zspacing];
    dim1 = [dim1 round(b(1)/config.div*i)];
    L1 = [L1 dim1(end)*config.pixel];
end
dim2 = [dim2 b(2)];
L2 = [L2 dim2(end)*config.zspacing];
dim1 = [dim1 b(1)];
L1 = [L1 dim1(end)*config.pixel];
% account for cropping
% crop_ind = [minrow maxrow mincol maxcol]
L2 = L2 + crop_index(3) * config.zspacing;
L1 = L1 + crop_index(1) * config.pixel;
%
L2c = {};
for i=1:length(L2)
    L2c = [L2c num2str(round(L2(i)),'%.0f')];
end
L1c = {};
for i=1:length(L1)
    L1c = [L1c num2str(round(L1(i)),'%.0f')];
end
ax=gca;
ax.XTick = dim2;
ax.XTickLabel = L2c;
ax.YTick = dim1;
ax.YTickLabel = L1c;
xlabel(xlbl);
ylabel(ylbl);
dim = get(gca, 'Position');
title({ [config.fname] },'fontsize',6,'interpreter','none');
end