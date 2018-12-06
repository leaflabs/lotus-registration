function bead_psf = axial_psf (config)

timestamp = datestr(datetime('now'));
fdiary = sprintf('%s%s.log',config.outpath,timestamp);
%fdiary = 'axial_psf.log'
diary(fdiary)
disp('%%');
disp(['%% ' datestr(datetime)]);
disp('%%');

config.pixel = 0.5; % um
config.zslices = 60;
config.zrange = [-300, 300];
% JPK
config.crop_y_pixel = ceil(config.crop_um/config.pixel);
config.crop_z_pixel = ceil(config.crop_um/config.zspacing);
config.crop_z = config.crop_z_pixel;
config.crop_row = config.crop_y_pixel;
config.div = 4;
%JPK
config.latBuf = 1;
config.axialBuf = 1;
config.threshold = 0.5;
config.figpos = [230        -895        1440         823];

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
    if ~isempty(find(config.blacklist==i))
        msg = sprintf('Skipping index %d since blacklisted.',i);
        disp(msg);
        continue;
    end
    %%% find z value of bead
    row = xycoords(i,2);
    col = xycoords(i,1);
    z = [z get_z( row, col, XguessSAVE1 )];
    if z(end) < 0
        disp('Error! No z value found.');
        continue
    end
    disp(sprintf('\n\n'));
    msg = sprintf('bead %d\ncentroid pixels [%d %d %d]',length(index),row,col,z(end));
    disp(msg);
    %%% slice volume at dim2=col to create a dim13 slice 
    slice13 = squeeze(XguessSAVE1(:,col,:));
    % crop_index = [minrow maxrow mincol maxcol]
    crop_index13 = cropIndex(config, size(slice13), row, z(end));
    if crop_index13<0
        continue;
    end
    cropped_slice13 = cropImag(slice13, crop_index13);
    %%% slice volume at dim3=z to create a dim12 slice 
    slice12 = squeeze(XguessSAVE1(:,:,z(end)));
    % crop_index = [minrow maxrow mincol maxcol]
    crop_index12 = cropIndex(config, size(slice12), row, col);
    if crop_index12<0
        continue;
    end
    cropped_slice12 = cropImag(slice12, crop_index12);
    %%% plot image
    f = figure;
    h = subplot(1,3,1);
    [dim, ztick_values, ztick_labels] = plot_slice_cropped ( config, h, ...
        cropped_slice13, size(slice13), crop_index13, ...
        ' z, i.e. depth (um)', ' y, i.e. height (um)');
    h = subplot(1,3,2);
    lateralVals = plot_psf_cropped(config, h, cropped_slice13, 1, 'y (um)', config.latBuf);
    h = subplot(1,3,3);
    axialVals = plot_psf_cropped(config, h, cropped_slice13, 2, 'z (um)', config.axialBuf);
    %
    keyboard
    % calculate PSFs
    latC = latPSF (config, f, lateralVals);
    axialC = axialPSF (config, f, axialVals, crop_index,...
        dim,size(cropped_slice13),z(end),0);
    %
    if latC(2)>0 & axialC(2)>0
        bead_psf = [bead_psf;[latC(1) latC(2) axialC(1) axialC(2)]];
        index = [index i];
    end
    figure(f);
    text(50.0,0.5,sprintf('bead %d',i),'FontSize',12);
    %
    str=sprintf('%s%s_bead%03d.png',config.outpath,config.fname(1:end-4),i);
    print(f,str,'-dpng');
    keyboard
end

f = figure();
plotAllPSFs(f, bead_psf, config.fpath);
str=sprintf('%s%s_all_beads.png',config.outpath,config.fname(1:end-4));
print(f,str,'-dpng');
disp(sprintf('# Output file = %s',str));

diary off

end


