function bead_psf = axial_psf (config)
load([config.fpath config.fname]);
load([config.xypath config.xyname]);
xycoords = round(xycoords);
a = size(xycoords);
if length(xycoords)==0
    diary off;
    return;
end
% for each bead
index = [];
z = [];
bead_psf = [];
for i=1:a(1)
    %i
%     if ~isempty(find(config.blacklist==i))
%         msg = sprintf('Skipping index %d since blacklisted.',i);
%         disp(msg);
%         continue;
%     end
    %%% find z value of bead
    row = xycoords(i,2);
    col = xycoords(i,1);
    z = [z get_z( row, col, XguessSAVE1 )];
    if z(end) < 0
        disp('Error! No z value found.');
        continue
    end
    disp(sprintf('\n\n'));
    msg = sprintf('bead %d\ncentroid pixels [%d %d %d]', ...
        i,row,col,z(end));
    disp(msg);
    
    %%% slice volume at dim2=col to create a dim13 slice 
    slice13 = squeeze(XguessSAVE1(:,col,:));
    % crop_index = [minrow maxrow mincol maxcol]
    crop_index13 = cropIndex(config, size(slice13), ...
        row, config.crop_row, ...
        z(end), config.crop_z);
    if crop_index13<0
        msg = sprintf('Skipping bead.');
        disp(msg);
        bead_psf = [bead_psf;[i -1 -1 -1 -1 -1 -1 -1 ]];
        continue;
    end
    cropped_slice13 = cropImag(slice13, crop_index13);
    
    %%% slice volume at dim3=z to create a dim12 slice 
    slice12 = squeeze(XguessSAVE1(:,:,z(end)));
    % crop_index = [minrow maxrow mincol maxcol]
    crop_index12 = cropIndex(config, size(slice12), ...
        row, config.crop_row, ...
        col, config.crop_row);
    if crop_index12<0
        msg = sprintf('Skipping bead.');
        disp(msg);
        bead_psf = [bead_psf;[i -1 -1 -1 -1 -1 -1 -1 ]];
        continue;
    end
    cropped_slice12 = cropImag(slice12, crop_index12);
    
    %%% plot image
    f = figure;
    set(f,'Position',config.figpos);
    
    h = subplot(2,3,1);
    plot_slice_cropped ( config, h, ...
        cropped_slice13, size(slice13), crop_index13, ...
        ' dim 3 or z (um)', ' dim 1 or y (um)', ...
        config.zspacing, config.pixel);
    h = subplot(2,3,2);
    fwhm1a = plot_psf_cropped(h, ...
        cropped_slice13, 2, 'dim 1 or y (um)', config.pixel);
    h = subplot(2,3,3);
    fwhm3 = plot_psf_cropped(h, ...
        cropped_slice13, 1, 'dim 3 or z (um)', config.zspacing);
    %
    h = subplot(2,3,4);
    plot_slice_cropped ( config, h, ...
        cropped_slice12, size(slice12), crop_index12, ...
        ' dim 2 or x (um)', ' dim 1 or y (um)', ...
        config.pixel, config.pixel);
    h = subplot(2,3,5);
    fwhm1b = plot_psf_cropped(h, ...
        cropped_slice12, 2, 'dim 1 or y (um)', config.pixel);
    h = subplot(2,3,6);
    fwhm2 = plot_psf_cropped(h, ...
        cropped_slice12, 1, 'dim 2 or x (um)', config.pixel);
    %
    bead_psf = [bead_psf;[i ...
        row*config.pixel col*config.pixel z(end)*config.zspacing ...
        fwhm1a fwhm1b fwhm2 fwhm3]];
    %
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1);
    str=sprintf('%s bead%03d',config.label,i);
    text(0.05,0.3,str,'FontSize',40,'Color',[0.8 0.8 0.8],'Rotation',90,'Interpreter','none');
    %
    str=sprintf('%s%s%s_bead%03d.png',config.outpath,config.label,config.fname(1:end-4),i);
    print(f,str,'-dpng');
    str=sprintf('%s%s%s_bead%03d.eps',config.outpath,config.label,config.fname(1:end-4),i);
    print(f,str,'-depsc');
    close all;
end
end


