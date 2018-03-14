%function bead_psf = axial_psf (fname, beadname, zspacing, crop_um, scale, blacklist, config)
function [bead_psf_DDLFM, bead_psf_LFM2] = axial_psf (s, blacklist, config)

% get bead coordinates
xyname = [s.folder '/' s.name]
load(xyname);
xycoords = round(xycoords);
a = size(xycoords);
if length(xycoords)==0
    diary off;
    return;
end

% get DDLFM volume
volname = [xyname(1:end-12) '.mat']
config.volname = volname;
load(config.volname);
DDLFM = XguessSAVE1;


% get LFM2 volume
if isempty(strfind(config.volname, '__'))
    disp('WTF?');
    keyboard
end
v = strfind(config.volname, '/');
ppath = config.volname(1:v(end-1));
psfpath = [ppath 'psf'];
if ~exist(psfpath,'dir')
    status = mkdir(psfpath);
    if status == 1
        disp(['Created folder: ' psfpath]);
    else
        disp(['Error attempting to create folder:' psfpath]);
        status
        exit;
    end
end
prefix = [ppath 'vertical/Reconstructed/'];
suffix = config.volname(v(end)+1:end);
v = strfind(suffix, '__');
lmf2_file = [prefix suffix(1:v-1) '.mat']
load(lmf2_file);
LFM2 = XguessSAVE1;

% for each bead
index = [];
z = [];
bead_psf_DDLFM = [];
bead_psf_LFM2 = [];
%for i=1:a(1)
for i=1:4

    if ~isempty(find(blacklist==i))
        msg = sprintf('Skipping index %d since blacklisted.',i);
        disp(msg);
        continue;
    end
    
    row = xycoords(i,2);
    col = xycoords(i,1);
    
    %% analyze DDLFM
    % change config
    crop_um = 20;
    config.zspacing = 0.5;
    crop_y_pixel = ceil(crop_um/config.pixel);
    crop_z_pixel = ceil(crop_um/config.zspacing);
    config.crop_z = crop_z_pixel;
    config.crop_row = crop_y_pixel;
    config.bit = ceil(log2(single(max(max(max(DDLFM))))));
    
    z = get_z( row, col, DDLFM );
    %z = z(end);
    if z< 0
        disp('Error! No z value found.');
        continue
    end
    disp(sprintf('\n\n'));
    msg = sprintf('DDLFM: bead %d found at row = %d, col = %d, z = %d',i,row,col,z);
    disp(msg);
    
    [bead_center_y_1, fwhm_y_1, bead_center_z, fwhm_z] ...
        = sliceYZ(DDLFM, row, col, z, config, sprintf('%s/%s_bead%03d_YZ.png',psfpath,suffix(1:end-4),i));
    [bead_center_y_2, fwhm_y_2, bead_center_x, fwhm_x] ...
        = sliceYX(DDLFM, row, col, z, config, sprintf('%s/%s_bead%03d_YX.png',psfpath,suffix(1:end-4),i));
    if fwhm_y_1 > 0 & fwhm_y_2 > 0 & fwhm_x > 0 & fwhm_z > 0
        bead_psf_DDLFM = [bead_psf_DDLFM;[z, bead_center_y_1, fwhm_y_1, bead_center_z, fwhm_z, bead_center_y_2, fwhm_y_2, bead_center_x, fwhm_x]];
        %index = [index i];
    end
    
    %% analyze LFM2
    % change config
    crop_um = 60;
    config.zspacing = 4.0;
    crop_y_pixel = ceil(crop_um/config.pixel);
    crop_z_pixel = ceil(crop_um/config.zspacing);
    config.crop_z = crop_z_pixel;
    config.crop_row = crop_y_pixel;
    config.bit = ceil(log2(single(max(max(max(LFM2))))));

    row = row + config.clip(1);
    col = col + config.clip(3);
    z = get_z( row, col, LFM2 );
    if z < 0
        disp('Error! No z value found.');
        continue
    end
    disp(sprintf('\n\n'));
    msg = sprintf('LFM2: looking for bead %d at row = %d, col = %d, z = %d',i,row,col,z);
    disp(msg);
    [bead_center_y_1, fwhm_y_1, bead_center_z, fwhm_z] ...
        = sliceYZ(LFM2, row, col, z, config, sprintf('%s/psf/%s_bead%03d_YZ.png',ppath,suffix(1:v-1),i));
 
    crop_um = 20;
    crop_y_pixel = ceil(crop_um/config.pixel);
    crop_z_pixel = ceil(crop_um/config.zspacing);
    config.crop_z = crop_z_pixel;
    config.crop_row = crop_y_pixel;
    
    [bead_center_y_2, fwhm_y_2, bead_center_x, fwhm_x] ...
        = sliceYX(LFM2, row, col, z, config, sprintf('%s/psf/%s_bead%03d_YX.png',ppath,suffix(1:v-1),i));
    if fwhm_y_1 > 0 & fwhm_y_2 > 0 & fwhm_x > 0 & fwhm_z > 0
        bead_psf_LFM2 = [bead_psf_LFM2;[z, bead_center_y_1, fwhm_y_1, bead_center_z, fwhm_z, bead_center_y_2, fwhm_y_2, bead_center_x, fwhm_x]];
    end
    close all;
end
end



function [bead_center_y, fwhm_y, bead_center_x, fwhm_x] = sliceYX(XguessSAVE1, row, col, z, config, beadname)
% slice volume
csec = squeeze(XguessSAVE1(:,:,z));
%
latBuf = 1;
%axialBuf = 1;
crop_index = cropIndexXY(config, size(csec), row, col);
if crop_index<0
    bead_center_y = -1;
    fwhm_y = -1;
    bead_center_x = -1;
    fwhm_x = -1;
    return;
end
csec_cropped = cropImag(csec, crop_index);
% plot image
h = figure;
[dim, ztick_values, ztick_labels] = plot_slice_croppedXY ( config, h, ...
    csec_cropped, size(csec), crop_index, beadname);
lateralVals = plot_lat_psf_cropped(config, h, csec_cropped);
axialVals = plot_axial_psf_cropped(config, h, csec_cropped);
% calculate PSFs
[bead_center_y fwhm_y] = latPSF (config, h, lateralVals);
%latC = latPSF (config, h, lateralVals);
[bead_center_x fwhm_x] = latPSFx (config, h, axialVals, crop_index,...
    dim,size(csec_cropped),z,0);
%axialC = axialPSF (config, h, axialVals, crop_index,...
%    dim,size(csec_cropped),z,0);
figure(h);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.05,0.97,beadname,'FontSize',10,'Color',[0 0 0],'Interpreter','none');
print(h,beadname,'-dpng');
end

function [bead_center_y, fwhm_y, bead_center_z, fwhm_z] = sliceYZ(XguessSAVE1, row, col, z, config, beadname)
% slice volume
csec = squeeze(XguessSAVE1(:,col,:));
%
latBuf = 1;
axialBuf = 1;
crop_index = cropIndex(config, size(csec), row, z);
if crop_index<0
    bead_center_y = -1;
    fwhm_y = -1;
    bead_center_z = -1;
    fwhm_z = -1;
    return;
end
csec_cropped = cropImag(csec, crop_index);
% plot image
h = figure;
[dim, ztick_values, ztick_labels] = plot_slice_cropped ( config, h, ...
    csec_cropped, size(csec), crop_index, beadname);
lateralVals = plot_lat_psf_cropped(config, h, csec_cropped);
axialVals = plot_axial_psf_cropped(config, h, csec_cropped);
% calculate PSFs
[bead_center_y fwhm_y] = latPSF (config, h, lateralVals);
%latC = latPSF (config, h, lateralVals);
[bead_center_z fwhm_z] = axialPSF (config, h, axialVals, crop_index,...
    dim,size(csec_cropped),z,0);
%axialC = axialPSF (config, h, axialVals, crop_index,...
%    dim,size(csec_cropped),z,0);
figure(h);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.05,0.97,beadname,'FontSize',10,'Color',[0 0 0],'Interpreter','none');
print(h,beadname,'-dpng');
end

function [bead_center, fwhm] = latPSFx (config, handle, axialVals, crop_index, dim, b, z, myflag)
myStr = {};
threshold = config.threshold;
N = length(axialVals)-1;
axialPosInt0 = 0;
axialPosInt1 = 0;
riseFlag = false;
fallFlag = false;
zvals = [z-5:z+5];
offset = crop_index(3) * config.pixel;
% for each axialVals = number of columns in cropped image
for j=1:length(axialVals)-1
    if j<config.axialBuf | N-j+1<config.axialBuf
        j
        continue;
    end
    if axialVals(j) <= threshold && axialVals(j+1) > threshold
        if myflag
            keyboard
        end
        if ~riseFlag
            riseFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max axial spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            axialPosInt0 = -1;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        %
        axialPos0 = offset + j*config.pixel;
        axialPos1 = axialPos0 + config.pixel;
        scale = ( threshold - axialVals(j) ) / ( axialVals(j+1) - axialVals(j) );
        axialPosInt0 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(axialPosInt0))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    elseif axialVals(j) >= threshold && axialVals(j+1) < threshold
        if myflag
            keyboard
        end
        if ~fallFlag
            fallFlag = true;
        else
            thisStr = ['Warning! Too many falling threshold crossings! Full-width half-max axial spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            axialPosInt0 = -1;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        axialPos0 = offset + j*config.pixel;
        axialPos1 = axialPos0 + config.pixel;
        scale = ( axialVals(j) - threshold) / ( axialVals(j) - axialVals(j+1) );
        axialPosInt1 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(axialPosInt1))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    end
end

if myflag
    keyboard
end
if axialPosInt0 < 0
    thisStr = ['Warning! Full-width half-max axial spread not calculated.'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    myStr = myStr(3:end);
    % calculate bead center
    bead_center = -1;
    fwhm = -1;
else
    if myflag
        keyboard
    end
    fwhm = axialPosInt1-axialPosInt0;
    thisStr = ['Full-width half-max axial spread = ' num2str(round(fwhm)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = 0.5*(axialPosInt1+axialPosInt0);
end
figure(handle);
% x = xlim;
% y = ylim;
% hold on;
% plot(xlim,(y(1)+config.latBuf-1)*[1 1],'r-');
% plot(xlim,(y(2)-config.latBuf+1)*[1 1],'r-');
% hold off;
% subplot(2,2,3);
% text(5.0,0.5,myStr,'FontSize',6);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.5,0.2,myStr,'FontSize',10,'Color',[0 0 0],'Interpreter','none');
%% duplicate
% thisdim = get(gca,'Position');
% set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
% set(gca,'XTickLabel',{});
% xlim([1 b(2)]);
%
%
% x = xlim;
% y = ylim;
% hold on;
% plot((x(1)+config.axialBuf-1)*[1 1],ylim,'r-');
% plot((x(2)-config.axialBuf+1)*[1 1],ylim,'r-');
% hold off;
end



function [bead_center, fwhm] = axialPSF (config, handle, axialVals, crop_index, dim, b, z, myflag)
myStr = {};
threshold = config.threshold;
N = length(axialVals)-1;
axialPosInt0 = 0;
axialPosInt1 = 0;
riseFlag = false;
fallFlag = false;
zvals = [z-5:z+5];
offset = crop_index(3) * config.zspacing;
% for each axialVals = number of columns in cropped image
for j=1:length(axialVals)-1
    if j<config.axialBuf | N-j+1<config.axialBuf
        j
        continue;
    end
    if axialVals(j) <= threshold && axialVals(j+1) > threshold
        if myflag
            keyboard
        end
        if ~riseFlag
            riseFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max axial spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            axialPosInt0 = -1;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        %
        axialPos0 = offset + j*config.zspacing;
        axialPos1 = axialPos0 + config.zspacing;
        scale = ( threshold - axialVals(j) ) / ( axialVals(j+1) - axialVals(j) );
        axialPosInt0 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(axialPosInt0))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    elseif axialVals(j) >= threshold && axialVals(j+1) < threshold
        if myflag
            keyboard
        end
        if ~fallFlag
            fallFlag = true;
        else
            thisStr = ['Warning! Too many falling threshold crossings! Full-width half-max axial spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            axialPosInt0 = -1;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        axialPos0 = offset + j*config.zspacing;
        axialPos1 = axialPos0 + config.zspacing;
        scale = ( axialVals(j) - threshold) / ( axialVals(j) - axialVals(j+1) );
        axialPosInt1 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(axialPosInt1))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    end
end

if myflag
    keyboard
end
if axialPosInt0 < 0
    thisStr = ['Warning! Full-width half-max axial spread not calculated.'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    myStr = myStr(3:end);
    % calculate bead center
    bead_center = -1;
    fwhm = -1;
else
    if myflag
        keyboard
    end
    fwhm = axialPosInt1-axialPosInt0;
    thisStr = ['Full-width half-max axial spread = ' num2str(round(fwhm)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = 0.5*(axialPosInt1+axialPosInt0);
end
figure(handle);
% x = xlim;
% y = ylim;
% hold on;
% plot(xlim,(y(1)+config.latBuf-1)*[1 1],'r-');
% plot(xlim,(y(2)-config.latBuf+1)*[1 1],'r-');
% hold off;
%subplot(2,2,3);
%text(5.0,0.5,myStr,'FontSize',6);

%% duplicate
% thisdim = get(gca,'Position');
% set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
% set(gca,'XTickLabel',{});
% xlim([1 b(2)]);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.5,0.2,myStr,'FontSize',10,'Color',[0 0 0],'Interpreter','none');
%
%
% x = xlim;
% y = ylim;
% hold on;
% plot((x(1)+config.axialBuf-1)*[1 1],ylim,'r-');
% plot((x(2)-config.axialBuf+1)*[1 1],ylim,'r-');
% hold off;
end

function [bead_center, fwhm] = latPSF (config, handle, lateralVals)
myStr = {};
threshold = config.threshold;
N = length(lateralVals)-1;
latPosInt0 = 0;
latPosInt1 = 0;
riseFlag = false;
fallFlag = false;
for j=1:length(lateralVals)-1
    if j<config.latBuf | N-j+1<config.latBuf
        j
        continue;
    end
    if lateralVals(j) <= threshold && lateralVals(j+1) > threshold
        if ~riseFlag
            riseFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max lat spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            latPosInt0 = -1;
        end
        if latPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        latPos0 =  j*config.pixel;
        latPos1 = (j+1)*config.pixel;
        scale = ( threshold - lateralVals(j) ) / ( lateralVals(j+1) - lateralVals(j) );
        latPosInt0 = (j+scale)*config.pixel;
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(latPosInt0))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    elseif lateralVals(j) >= threshold && lateralVals(j+1) < threshold
        if ~fallFlag
            fallFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max lat spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            latPosInt0 = -1;
        end
        if latPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        latPos0 =  j*config.pixel;
        latPos1 = (j+1)*config.pixel;
        scale = ( lateralVals(j) - threshold) / ( lateralVals(j) - lateralVals(j+1) );
        latPosInt1 = (j+scale)*config.pixel;
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(latPosInt1))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    end
end
if latPosInt0 < 0
    thisStr = ['Warning! Full-width half-max lat spread not calculated.'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    myStr = myStr(3:end);
    % calculate bead center
    bead_center = -1
    fwhm = -1;
else
    fwhm = latPosInt1-latPosInt0;
    thisStr = ['Full-width half-max lateral spread = ' num2str(round(fwhm)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = 0.5*(latPosInt1+latPosInt0);
end

figure(handle);
%subplot(2,2,2);
%text(0.1,15.0,myStr,'FontSize',6);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.55,0.45,myStr,'FontSize',10,'Color',[0 0 0],'Interpreter','none');
end

function crop_ind = cropIndexXY(config, a, row, col)
dy = ceil(config.crop_row/2);
dx = dy;
ylim = [row-dy row+dy];
xlim = [col-dx col+dx];
ylim(1) = boundsCheck (ylim(1)<1, ylim(1), 1, 'Min row', a);
ylim(2) = boundsCheck (ylim(2)>a(1), ylim(2), a(1), 'Max row', a);
xlim(1) = boundsCheck (xlim(1)<1, xlim(1), 1, 'Min col', a);
xlim(2) = boundsCheck (xlim(2)>a(2), xlim(2), a(2), 'Max col', a);
if find([ylim(1) ylim(2) xlim(1) xlim(2)] < 0)
    crop_ind = -1;
    return
end
crop_ind = [ylim(1) ylim(2) xlim(1) xlim(2)];
end


function crop_ind = cropIndex(config, a, row, z)
dr = ceil(config.crop_row/2);
dz = ceil(config.crop_z/2);
rlim = [row-dr row+dr];
zlim = [z-dz z+dz];
rlim(1) = boundsCheck (rlim(1)<1, rlim(1), 1, 'Min row', a);
rlim(2) = boundsCheck (rlim(2)>a(1), rlim(2), a(1), 'Max row', a);
zlim(1) = boundsCheck (zlim(1)<1, zlim(1), 1, 'Min z', a);
zlim(2) = boundsCheck (zlim(2)>a(2), zlim(2), a(2), 'Max z', a);
if find([rlim(1) rlim(2) zlim(1) zlim(2)] < 0)
    crop_ind = -1;
    return
end
crop_ind = [rlim(1) rlim(2) zlim(1) zlim(2)];
end


function csec_cropped = cropImag(csec, crop_index)
csec_cropped = csec(crop_index(1):crop_index(2),crop_index(3):crop_index(4));
a = size(csec);
msg = sprintf('\nOrig slice is %d row by %d z',a(1),a(2));
disp(msg);
b = size(csec_cropped);
msg = sprintf('Cropped slice is %d row by %d z',b(1),b(2));
disp(msg);
end

function out = boundsCheck (cond, val, bound, mystr, a)
if cond
    msg = sprintf('Warning. %s pixel %d exceeds bound %d and is outside slice (%d,%d).', ...
        mystr,val, bound, a(1),a(2));
    disp(msg);
    msg = sprintf('Skipping bead.');
    disp(msg);
    out = -1;
else
    out = val;
end
end


function axialVals = plot_axial_psf_cropped (config, handle, csec)
subplot(2,2,3);
a = max(csec);
axialVals = single(a)/max(single(a(config.axialBuf:length(a)-config.axialBuf+1)));
plot(axialVals,'*-');
set(gca, 'YDir', 'reverse');
set(gca,'XTickLabel',{});
b = size(csec);
xlim([1 b(2)]);
ylabel('Normalized intensity');
ylim([0 1]);
end

function lateralVals = plot_lat_psf_cropped ( config, handle, csec )
a = max(csec');
lateralVals = single(a)/max(single(a(config.latBuf:length(a)-config.latBuf+1)));
figure(handle);
subplot(2,2,2);
plot(lateralVals,[1:length(a)],'*-');
set(gca, 'YDir', 'reverse');
set(gca,'YTickLabel',{});
xlabel('Normalized intensity');
ylim([0 length(a)]);
xlim([0 1]);
end


function [dim,z,Lzc] = plot_slice_croppedXY ( config, handle, ...
    csec_cropped, a, crop_index, beadname )
figure(handle);
h = subplot(2,2,1);
im = imagesc(csec_cropped,[0 2^config.bit]);
%
% transform figure so that plot axes have same scale
% to do so,
% the target crop size is config.crop_um, e.g. 100 um
% given crop_index calculate actual size in um
% get size of plot (in whatever units)
% scale plot so 1 um in x is same size on screen as 1 um in y
%
b = size(csec_cropped);
pixel = config.pixel;
vertical_microns = single(b(1)) * pixel;
depth_microns = (single(b(2))-1.0) * single(config.pixel);
disp(['image size = ' num2str(vertical_microns) ' um tall x ' num2str(depth_microns) ' um wide' ]);
% here, width is depth, and height is vertical
pos = get(h,'Position'); %  gives x left, y bottom, width, height
newpos = pos;
% scale width
newpos(3) = pos(3) * vertical_microns / depth_microns;
set(h,'Position', newpos);
%
% relabel axes to microns
% to do so,
% specify XTick = tick locations in csec_cropped space
% specify XTickLabel = labels for each tick
%
% relative to cropped image (i.e. origin in upper left corner)
z = [1];
Lz = [1];
y = [1];
Ly = [1];
for i=1:config.div-1
    z = [z round(b(2)/config.div*i)];
    Lz = [Lz z(end)*config.pixel];
    y = [y round(b(1)/config.div*i)];
    Ly = [Ly y(end)*config.pixel];
end
z = [z b(2)];
Lz = [Lz z(end)*config.pixel];
y = [y b(1)];
Ly = [Ly y(end)*config.pixel];
% account for cropping
Lz = Lz + crop_index(3) * config.pixel;
Ly = Ly + crop_index(1) * config.pixel;
Lzc = {};
for i=1:length(Lz)
    Lzc = [Lzc num2str(round(Lz(i)),'%.0f')];
end
Lyc = {};
for i=1:length(Ly)
    Lyc = [Lyc num2str(round(Ly(i)),'%.0f')];
end

ax=gca;
ax.XTick = z;
ax.XTickLabel = Lzc;
ax.YTick = y;
ax.YTickLabel = Lyc;

xlabel(' x or width (um)');
ylabel(' y or height (um)');
dim = get(gca, 'Position');
%title(beadname,'fontsize',6,'interpreter','none');
end


function [dim,z,Lzc] = plot_slice_cropped ( config, handle, ...
    csec_cropped, a, crop_index, beadname )
figure(handle);
h = subplot(2,2,1);
im = imagesc(csec_cropped,[0 2^config.bit]);
%
% transform figure so that plot axes have same scale
% to do so,
% the target crop size is config.crop_um, e.g. 100 um
% given crop_index calculate actual size in um
% get size of plot (in whatever units)
% scale plot so 1 um in x is same size on screen as 1 um in y
%
b = size(csec_cropped);
pixel = config.pixel;
vertical_microns = single(b(1)) * pixel;
depth_microns = (single(b(2))-1.0) * single(config.zspacing);
disp(['image size = ' num2str(vertical_microns) ' um tall x ' num2str(depth_microns) ' um wide' ]);
% here, width is depth, and height is vertical
pos = get(h,'Position'); %  gives x left, y bottom, width, height
newpos = pos;
% scale width
newpos(3) = newpos(4) * depth_microns / vertical_microns;
set(h,'Position', newpos);
%
% relabel axes to microns
% to do so,
% specify XTick = tick locations in csec_cropped space
% specify XTickLabel = labels for each tick
%
% relative to cropped image (i.e. origin in upper left corner)
z = [1];
Lz = [1];
y = [1];
Ly = [1];
for i=1:config.div-1
    z = [z round(b(2)/config.div*i)];
    Lz = [Lz z(end)*config.zspacing];
    y = [y round(b(1)/config.div*i)];
    Ly = [Ly y(end)*config.pixel];
end
z = [z b(2)];
Lz = [Lz z(end)*config.zspacing];
y = [y b(1)];
Ly = [Ly y(end)*config.pixel];
% account for cropping
Lz = Lz + crop_index(3) * config.zspacing;
Ly = Ly + crop_index(1) * config.pixel;
Lzc = {};
for i=1:length(Lz)
    Lzc = [Lzc num2str(round(Lz(i)),'%.0f')];
end
Lyc = {};
for i=1:length(Ly)
    Lyc = [Lyc num2str(round(Ly(i)),'%.0f')];
end

ax=gca;
ax.XTick = z;
ax.XTickLabel = Lzc;
ax.YTick = y;
ax.YTickLabel = Lyc;

xlabel(' z or depth (um)');
ylabel(' y or height (um)');
dim = get(gca, 'Position');
%title(beadname,'fontsize',6,'interpreter','none');
end


function z = get_z ( row, col, mymatrix )
a = size(mymatrix);
if row > a(1)
    msg = sprintf('Skipping bead, row %d is outside of image (%d)',row,a(1));
    disp(msg);
    z=-1;
    return;
elseif col > a(2)
    msg = sprintf('Skipping bead, col %d is outside of image (%d)',col,a(2));
    disp(msg);
    z=-1;
    return;
end
mymax = -1;
z = -1;
for i=1:a(3)
    if row<1 | col<1 | i<1
        disp('WTF?');
        keyboard;
    end
    val = mymatrix(row,col,i);
    if val > mymax
        mymax = val;
        z = i;
    end
end
end
