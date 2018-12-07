function plot_slice_cropped ( config, handle, ...
    csec_cropped, a, crop_index, xlbl, ylbl, scalex, scaley)
%%% render slice as image
subplot(handle);
bits = ceil(log2(single(max(max(csec_cropped)))));
cmap = gray(2^bits);
colormap(cmap);
imagesc(csec_cropped,[0 2^bits]);
daspect([scaley, scalex, 1]);
xlabel(xlbl);
ylabel(ylbl);
%
% crop_index = [minrow maxrow mincol maxcol]
labels_numX = linspace(crop_index(3)*scalex, crop_index(4)*scalex, config.div+1);
labelsX = {};
for i=1:length(labels_numX)
    labelsX = [labelsX num2str(round(labels_numX(i)),'%.0f')];
end
x = xlim;
ticksX = linspace(x(1), x(2), config.div+1);
labels_numY = linspace(crop_index(1)*scaley, crop_index(2)*scaley, config.div+1);
labelsY = {};
for i=1:length(labels_numX)
    labelsY = [labelsY num2str(round(labels_numY(i)),'%.0f')];
end
y = ylim;
ticksY = linspace(y(1), y(2), config.div+1);
%
ax=gca;
ax.XTick = ticksX;
ax.XTickLabel = labelsX;
ax.YTick = ticksY;
ax.YTickLabel = labelsY;
ax.FontSize = 30;
end