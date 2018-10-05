function save_2d_max_projections_compact (LFM1, LFM2, new, param, flag, str)

f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1);
yz = squeeze(max(LFM1,[],2));
xy = squeeze(max(LFM1,[],3));
xz = squeeze(max(LFM1,[],1))';
a = size(LFM1);
big_image = [a(3)+a(1) a(3)+a(2)];
imagesc(zeros(big_image));
a = size(yz);
x = [1 a(2)];
y = [1 a(1)];
imagesc('XData',x,'YData',y,'CData',fliplr(yz));
b = size(xy);
x = [a(2)+1 b(2)+a(2)];
y = [1 b(1)];
imagesc('XData',x,'YData',y,'CData',xy);
c = size(xz);
x = [a(2)+1 b(2)+a(2)];
y = [b(1)+1 b(1)+c(1)];
imagesc('XData',x,'YData',y,'CData',xz);
title('LFM1');
colorbar();
xlabel('pixels');
ylabel('pixels');
daspect([1,1,1]);

subplot(1,3,2);
yz = squeeze(max(LFM2,[],2));
xy = squeeze(max(LFM2,[],3));
xz = squeeze(max(LFM2,[],1))';
a = size(LFM2);
big_image = [a(3)+a(1) a(3)+a(2)];
imagesc(zeros(big_image));
a = size(yz);
x = [1 a(2)];
y = [1 a(1)];
imagesc('XData',x,'YData',y,'CData',fliplr(yz));
b = size(xy);
x = [a(2)+1 b(2)+a(2)];
y = [1 b(1)];
imagesc('XData',x,'YData',y,'CData',xy);
c = size(xz);
x = [a(2)+1 b(2)+a(2)];
y = [b(1)+1 b(1)+c(1)];
imagesc('XData',x,'YData',y,'CData',xz);
title('LFM2');
colorbar();
xlabel('pixels');
ylabel('pixels');
daspect([1,1,1]);

if flag
    subplot(1,3,3);
    yz = squeeze(max(new,[],2));
    xy = squeeze(max(new,[],3));
    xz = squeeze(max(new,[],1))';
    a = size(new);
    big_image = [a(3)+a(1) a(3)+a(2)];
    imagesc(zeros(big_image));
    a = size(yz);
    x = [1 a(2)];
    y = [1 a(1)];
    imagesc('XData',x,'YData',y,'CData',fliplr(yz));
    b = size(xy);
    x = [a(2)+1 b(2)+a(2)];
    y = [1 b(1)];
    imagesc('XData',x,'YData',y,'CData',xy);
    c = size(xz);
    x = [a(2)+1 b(2)+a(2)];
    y = [b(1)+1 b(1)+c(1)];
    imagesc('XData',x,'YData',y,'CData',xz);
    title('DDLFM');
    colorbar();
    xlabel('pixels');
    ylabel('pixels');
    daspect([1,1,1]);
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.1,0.98,['LFM1 = ' param.inputFilePath1 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
text(0.1,0.96,['LFM2 = ' param.inputFilePath2 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
if flag
    text(0.1,0.94,['DDLFM = ' param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
end
text(0.4,0.1,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
save_plot(f, str);

end
