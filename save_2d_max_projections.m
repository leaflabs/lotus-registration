function save_2d_max_projections (LFM1, LFM2, new, param, flag, str)

f = figure('units','normalized','outerposition',[0 0 1 1]);
yz = squeeze(max(LFM1,[],2));
subplot(2,6,1);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM1');
%colorbar();
daspect([1,1,1]);
set(gca,'XDir','reverse');

xy = squeeze(max(LFM1,[],3));
subplot(2,6,2);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM1');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(LFM1,[],1))';
subplot(2,6,8);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
title('LFM1');
colorbar();
daspect([1,1,1]);

yz = squeeze(max(LFM2,[],2));
subplot(2,6,3);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM2');
%colorbar();
daspect([1,1,1]);
set(gca,'XDir','reverse');


xy = squeeze(max(LFM2,[],3));
subplot(2,6,4);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM2');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(LFM2,[],1))';
subplot(2,6,10);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
title('LFM2');
colorbar();
daspect([1,1,1]);

if flag
    yz = squeeze(max(new,[],2));
    subplot(2,6,5);
    dr = ceil(log2(single(max(max(yz)))));
    if dr == 0 || max(max(max(new)))==0
        fprintf('Warning! dr =0.\n');
        dr = 8;
    end
    imagesc(yz,[0 2^dr]);
    xlabel('three [pixels]');
    ylabel('one [pixels]');
    title('DDLFM');
    %colorbar();
    daspect([1,1,1]);
    set(gca,'XDir','reverse');
    
    xy = squeeze(max(new,[],3));
    subplot(2,6,6);
    dr = ceil(log2(single(max(max(xy)))));
    if dr == 0 || max(max(max(new)))==0
        fprintf('Warning! dr =0.\n');
        dr = 8;
    end
    imagesc(xy,[0 2^dr]);
    xlabel('two [pixels]');
    ylabel('one [pixels]');
    title('DDLFM');
    %colorbar();
    daspect([1,1,1]);
    
    xz = squeeze(max(new,[],1))';
    subplot(2,6,12);
    dr = ceil(log2(single(max(max(xz)))));
    if dr == 0 || max(max(max(new)))==0
        fprintf('Warning! dr =0.\n');
        dr = 8;
    end
    imagesc(xz,[0 2^dr]);
    xlabel('two [pixels]');
    ylabel('three [pixels]');
    title('DDLFM');
    colorbar();
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
