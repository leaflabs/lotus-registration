function getbeads (fvec, outpath)
for i=1:length(fvec)
    
    fname = [fvec(i).folder '/' fvec(i).name]
    load(fname);
    if exist('Xvolume','var')
        xy = squeeze(sum(Xvolume,3));
    elseif exist('XguessSAVE1','var')
        xy = squeeze(sum(XguessSAVE1,3));
    else
        disp('WTF?! Unknown data name.');
        keyboard
    end
    
    figure(1);
    dynamic_range = round(log2(max(max(xy))));
    imagesc(xy,[0 2^dynamic_range]);
    xlabel('x');
    ylabel('y');
    colorbar();
    xycoords = [];
    hold on;
    while (1)
        a = ginput(1)
        if a(1)<50 & a(2) <50
            break;
        end
        xycoords = [xycoords;a];
        plot(a(1),a(2),'k*');
    end
    hold off;
    
    fout = [outpath fvec(i).name(1:end-4) '_bead_xy.mat']
    save(fout,'xycoords');
    
    fout = [outpath fvec(i).name(1:end-4) '_bead_xy.png']
    ff=getframe(gcf);
    [X, map] = frame2im(ff);
    imwrite(X, fout);
end
end