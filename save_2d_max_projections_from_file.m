function save_2d_max_projections_from_file (mypath, filename)

fname = [mypath '/' filename];
fprintf('\nLoading volume %s\n\n',fname);
LFM = loadData2(fname);

f = figure('units','normalized','outerposition',[0 0 1 1]);
yz = squeeze(max(LFM,[],2));
subplot(2,2,1);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM2 registered');
%colorbar();
daspect([1,1,1]);
set(gca,'XDir','reverse');

xy = squeeze(max(LFM,[],3));
subplot(2,2,2);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM2 registered');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(LFM,[],1))';
subplot(2,2,4);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
title('LFM2 registered');
colorbar();
daspect([1,1,1]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.1,0.98,['LFM = ' fname],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
fout=sprintf('%s/%s_MIP.png',mypath,filename(1:end-4));
print(f,fout,'-dpng');

end

function out = loadData2 (f)
if exist(f,'file') == 2
    load(f);
    if exist('out','var')
        if isa(out,'uint8')
            return;
        else
            disp('Warning input data is not recognized.');
            keyboard
        end
    else
        disp('Warning input data is not recognized.');
        keyboard
    end
else
    fprintf('File not found:\n%s\n\n',f);
    keyboard;
end
end