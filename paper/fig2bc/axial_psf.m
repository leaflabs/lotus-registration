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
    %lateralVals = plot_lat_psf_cropped(config, h, cropped_slice13);
    lateralVals = plot_psf_cropped(config, h, cropped_slice13, 1, 'y (um)', config.latBuf);
    h = subplot(1,3,3);
    axialVals = plot_psf_cropped(config, h, cropped_slice13, 2, 'z (um)', config.axialBuf);
    %axialVals = plot_axial_psf_cropped(config, h, cropped_slice13);
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

function out = plotAllPSFs (handle, bead_center, fpath)
subplot(3,1,2);
% for i=1:a(1)
% x axis = z pos = bead_center(:,3)
% y axis = fwhm in z = bead_center(:,4)
% n=[7 6 1:4];
% plot(bead_center(n,3),bead_center(n,4),'-o');
plot(bead_center(:,3),bead_center(:,4),'o');
% %
xlabel('z, i.e. depth (um)');
ylabel('full-width half-max in z (um)');
% plot i's
% tweak = 1.00;
% %for i=n
text(0, 70,['N = ' num2str(length(bead_center)) ' beads'],'FontSize',10);
% %    text(tweak*double(bead_center(i,3)),tweak*double(bead_center(i,4)),sprintf('[%d,%d]',zrange(i,1),zrange(i,2)),'FontSize',10);
% %end
% %hold off;
subplot(3,1,3);
% % x axis = y pos = height = bead_center(:,1)
% % y axis = fwhm in y = bead_center(:,2)
% plot(bead_center(n,3),bead_center(n,2),'-o');
plot(bead_center(:,3),bead_center(:,2),'o');
% %
xlabel('z, i.e. depth (um)');
ylabel('full-width half-max in y (um)');
% %
% % plot i's
% %for i=n
% %    text(tweak*double(bead_center(i,3)),tweak*double(bead_center(i,2)),sprintf('[%d,%d]',zrange(i,1),zrange(i,2)),'FontSize',10);
% %end
subplot(3,1,1);
% % x axis = y pos = height = bead_center(:,1)
% % y axis = fwhm in y = bead_center(:,2)
% plot(bead_center(n,3),bead_center(n,2),'-o');
plot(bead_center(:,4),bead_center(:,2),'o');
% %
xlabel('full-width half-max in z (um)');
ylabel('full-width half-max in y (um)');
%axis equal;

% figure(handle);
% str=sprintf('%ssummary.png',fpath);
% f=getframe(gcf);
% [X, map] = frame2im(f);
% imwrite(X, str);
%disp(sprintf('# Output file = %s',str));
end



function bead_center = axialPSF (config, handle, axialVals, crop_index, dim, b, z, myflag)
myStr = {};
threshold = config.threshold;
N = length(axialVals)-1;
axialPosInt0 = 0;
axialPosInt1 = 0;
riseFlag = false;
fallFlag = false;
zvals = [z-5:z+5];
% JPK
offset = crop_index(3) * config.zspacing;
%
% for each axialVals = number of columns in cropped image
for j=1:length(axialVals)-1
    axialVals(j);
    if j<config.axialBuf | N-j+1<config.axialBuf
        %j
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
%             % calculate bead center
%             bead_center = [-1 -1];
%             thisStr = ['Setting bead axial center and full-width half-max axial spread to %f and %f.', ...
%                 bead_center(1),bead_center(2)];
%             disp(thisStr);
%             return;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        %
        % crop_ind = [minrow maxrow mincol maxcol]
        axialPos0 = offset + j*config.zspacing;
        %axialPos0 = allLabels(zvals(j));
        axialPos1 = axialPos0 + config.zspacing;
        scale = ( threshold - axialVals(j) ) / ( axialVals(j+1) - axialVals(j) );
        axialPosInt0 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        %axialPosInt0 = (j+scale)*pixel;
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
%             % calculate bead center
%             bead_center = [-1 -1];
%             thisStr = ['Setting bead axial center and full-width half-max axial spread to %f and %f.', ...
%                 bead_center(1),bead_center(2)];
%             disp(thisStr);
%             return;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        axialPos0 = offset + j*config.zspacing;
        %axialPos0 = allLabels(zvals(j));
        axialPos1 = axialPos0 + config.zspacing;
        %axialPos0 =  allLabels(zvals(j));
        %axialPos1 = allLabels(zvals(j+1));
        scale = ( axialVals(j) - threshold) / ( axialVals(j) - axialVals(j+1) );
        axialPosInt1 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        %axialPosInt1 = (j+scale)*pixel;
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
    bead_center = [-1 -1];
    %thisStr = ['Setting bead axial center and full-width half-max axial spread to %f and %f.', ...
    %    bead_center(1),bead_center(2)];
    %disp(thisStr);
else
    if myflag
            keyboard
        end
    fwhmA = axialPosInt1-axialPosInt0;
    thisStr = ['Full-width half-max axial spread = ' num2str(round(fwhmA)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = [0.5*(axialPosInt1+axialPosInt0) fwhmA];
end
figure(handle);
x = xlim;
y = ylim;
hold on;
plot(xlim,(y(1)+config.latBuf-1)*[1 1],'r-');
plot(xlim,(y(2)-config.latBuf+1)*[1 1],'r-');
hold off;
subplot(2,2,3);
text(5.0,0.5,myStr,'FontSize',6);
%% duplicate
thisdim = get(gca,'Position');
set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
set(gca,'XTickLabel',{});
xlim([1 b(2)]);
%
%
x = xlim;
y = ylim;
hold on;
plot((x(1)+config.axialBuf-1)*[1 1],ylim,'r-');
plot((x(2)-config.axialBuf+1)*[1 1],ylim,'r-');
hold off;
end

function bead_center = latPSF (config, handle, lateralVals)
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
            % calculate bead center
%             bead_center = [-1 -1];
%             thisStr = ['Setting bead lat center and full-width half-max lat spread to %f and %f.', ...
%                 bead_center(1),bead_center(2)];
%             disp(thisStr);
%             return;
        end
        if latPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        latPos0 =  j*config.pixel;
        latPos1 = (j+1)*config.pixel;
        scale = ( threshold - lateralVals(j) ) / ( lateralVals(j+1) - lateralVals(j) );
        %latPosInt = latPos0 + (latPos1-latPos0)*scale; Equivalently...
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
%             % calculate bead center
%             bead_center = [-1 -1];
%             thisStr = ['Setting bead lat center and full-width half-max lat spread to %f and %f.', ...
%                 bead_center(1),bead_center(2)];
%             disp(thisStr);
%             return;
        end
        if latPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        latPos0 =  j*config.pixel;
        latPos1 = (j+1)*config.pixel;
        scale = ( lateralVals(j) - threshold) / ( lateralVals(j) - lateralVals(j+1) );
        %latPosInt = latPos0 + (latPos1-latPos0)*scale; Equivalently...
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
    bead_center = [-1 -1];
%     thisStr = ['Setting bead lat center and full-width half-max lat spread to %f and %f.', ...
%         bead_center(1),bead_center(2)];
%     disp(thisStr);
else
    fwhmL = latPosInt1-latPosInt0;
    thisStr = ['Full-width half-max lateral spread = ' num2str(round(fwhmL)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = [0.5*(latPosInt1+latPosInt0) fwhmL];
end

figure(handle);
subplot(2,2,2);
text(0.1,15.0,myStr,'FontSize',6);
%

end

function crop_ind = cropIndex(config, a, row, z)
dr = ceil(config.crop_row/2);
dz = ceil(config.crop_z/2);
% dr = (config.crop_row - 1) / 2;
% dz = (config.crop_z - 1) / 2;
rlim = [row-dr row+dr];
zlim = [z-dz z+dz];
% check bounds
%a = size(csec);
% if rlim(1) < 1 || rlim(2) > a(1) ...
%         || zlim(1) < 1 || zlim(2) > a(2)
%         csec_cropped = -1;
%     return
% end
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
msg = sprintf('Orig slice is %d by %d',a(1),a(2));
disp(msg);
b = size(csec_cropped);
msg = sprintf('Cropped slice is %d by %d',b(1),b(2));
disp(msg);
end

function out = boundsCheck (cond, val, bound, mystr, a)
if cond
    msg = sprintf('Warning. %s pixel %d exceeds bound %d and is outside slice (%d,%d).', ...
        mystr,val, bound, a(1),a(2));
    disp(msg);
    %out = bound;
    %msg = sprintf('%s pixel set to %d.', mystr, out);
    msg = sprintf('Skipping bead.');
    disp(msg);
    out = -1;
else
    out = val;
end
end

% function out = plot_axial_psf (config, handle, csec, axialBuf, pixel )
% subplot(2,2,3);
% a = max(csec);
% axialVals = single(a)/max(single(a(axialBuf:length(a)-axialBuf+1)));
% plot(axialVals,'*-');
% set(gca, 'YDir', 'reverse');
% %thisdim = get(gca,'Position');
% %set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
% set(gca,'XTickLabel',{});
% xlim([1 2*config.zslices+1]);
% ylabel('Normalized intensity');
% ylim([0 1]);
% end


% function axialVals = plot_axial_psf_cropped (config, handle, csec)
% subplot(2,2,3);
% a = max(csec);
% axialVals = single(a)/max(single(a(config.axialBuf:length(a)-config.axialBuf+1)));
% plot(axialVals,'*-');
% set(gca, 'YDir', 'reverse');
% %thisdim = get(gca,'Position');
% %set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
% set(gca,'XTickLabel',{});
% b = size(csec);
% xlim([1 b(2)]);
% ylabel('Normalized intensity');
% ylim([0 1]);
% end


function psfVals = plot_psf_cropped (config, handle, csec, dim, str, Buf)
subplot(handle);
a = max(csec,[],dim);
psfVals = single(a)/max(single(a(Buf:length(a)-Buf+1)));
vec = [1:length(psfVals)]*config.pixel;
plot(vec,psfVals,'*-');
%set(gca, 'YDir', 'reverse');
%thisdim = get(gca,'Position');
%set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
%set(gca,'XTickLabel',{});
%b = size(csec);
%xlim([1 b(2)]);
ylabel('Normalized intensity');
ylim([0 1]);
xlabel(str);
end


function lateralVals = plot_lat_psf_cropped ( config, handle, csec )
a = max(csec');
lateralVals = single(a)/max(single(a(config.latBuf:length(a)-config.latBuf+1)));
figure(handle);
subplot(2,2,2);
plot(lateralVals,[1:length(a)],'*-');
set(gca, 'YDir', 'reverse');
%thisdim = get(gca,'Position');
%set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) dim(4)]);
%ax=gca;
%ax.YTickLabel = {round(ax.YTick*config.pixel)};
%ylabel(' y, i.e. height (um)');
set(gca,'YTickLabel',{});
xlabel('Normalized intensity');
ylim([0 length(a)]);
xlim([0 1]);
end


% function out = plot_lat_psf ( handle, csec, latBuf, pixel )
% figure(handle);
% subplot(2,2,2);
% a = max(csec');
% lateralVals = single(a)/max(single(a(latBuf:length(a)-latBuf+1)));
% plot(lateralVals,[1:length(a)],'*-');
% set(gca, 'YDir', 'reverse');
% %thisdim = get(gca,'Position');
% %set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) dim(4)]);
% ax=gca;
% ax.YTickLabel = {round(ax.YTick*pixel)};
% %ylabel(' y, i.e. height (um)');
% set(gca,'YTickLabel',{});
% xlabel('Normalized intensity');
% ylim([0 length(a)]);
% xlim([0 1]);
% end

% function pixel = plot_slice ( config, handle, csec, z, row )
% figure(handle);
% subplot(2,2,1);
% imagesc(csec,[0 2^16]);
% colorbar();
% %
% % mark bead
% hold on;
% plot(z,row,'ow');
% hold off;
% %
% % transform figure so that plot axes have same scale
% a = size(csec);
% pixel = config.microlens_array_pitch / single(config.Nnum) / config.mag;
% vertical_microns = single(a(1)) * pixel;
% depth_microns = (single(a(2))-1.0) * single(config.zspacing);
% disp(['image size = ' num2str(vertical_microns) ' um tall x ' num2str(depth_microns) ' um wide' ]);
% pos = get(gcf,'Position'); %  gives x left, y bottom, width, height
% newpos = pos;
% newpos(3) = pos(4) * depth_microns / vertical_microns;
% set(gcf,'Position', newpos);
% % relabel axes to microns
% ax=gca;
% nz = config.zslices;
% dz = config.zspacing;
% allLabels = [-nz*dz:dz:nz*dz];
% %allLabels = [config.zrange(1):dz:config.zrange(2)];
% step = 10;
% tmp = [step:step:nz];
% ax.XTick = [-fliplr(tmp) 0 tmp] + (nz+1)*ones(1,2*length(tmp)+1);
% ax.XTickLabel = {allLabels(ax.XTick)};
% xlabel(' z, i.e. depth (um)');
% ax.YTickLabel = {round(ax.YTick*pixel)};
% ylabel(' y, i.e. height (um)');
% dim = get(gca, 'Position');
% % add title
% title({ [config.PSFname] ; ['Nnum = ' num2str(config.Nnum) ', maxIter = ' num2str(config.maxIter)]},'fontsize',6);
% end

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

% function [pixel,allLabels,dim,nz] = plot_slice_cropped ( config, handle, csec_cropped, z, row )
% figure(handle);
% subplot(2,2,1);
% imagesc(csec_cropped,[0 2^8]);
% colorbar();
% %
% % mark bead
% hold on;
% plot(z,row,'ow');
% hold off;
% %
% % transform figure so that plot axes have same scale
% a = size(csec_cropped);
% %pixel = config.microlens_array_pitch / single(config.Nnum) / config.mag;
% pixel = config.pixel;
% vertical_microns = single(a(1)) * pixel;
% depth_microns = (single(a(2))-1.0) * single(config.zspacing);
% disp(['image size = ' num2str(vertical_microns) ' um tall x ' num2str(depth_microns) ' um wide' ]);
% pos = get(gcf,'Position'); %  gives x left, y bottom, width, height
% newpos = pos;
% newpos(3) = pos(4) * depth_microns / vertical_microns;
% keyboard
% set(gcf,'Position', newpos);
% keyboard
% % relabel axes to microns
% ax=gca;
% nz = config.zslices;
% dz = config.zspacing;
% allLabels = [-nz*dz:dz:nz*dz];
% %allLabels = [config.zrange(1):dz:config.zrange(2)];
% %step = 3;
% %tmp = [step:step:nz];
% %ax.XTick = [-fliplr(tmp) 0 tmp] + (nz+1)*ones(1,2*length(tmp)+1);
% % TODO FIX HARDCODING
% ax.XTick = [1:11];
% labelindex = [z-5:z+5];
% % if isempty(find(labelindex>0))
% %     disp('WTF?')
% %     keyboard
% % end
% % if z==65
% %     keyboard
% % end
% keyboard
% labelindex = labelindex(find(labelindex>0));
% ax.XTickLabel = {allLabels(labelindex)};
% xlabel(' z, i.e. depth (um)');
% ax.YTickLabel = {round(ax.YTick*pixel)};
% ylabel(' y, i.e. height (um)');
% dim = get(gca, 'Position');
% % add title
% title({ [config.fname] },'fontsize',6,'interpreter','none');
% end

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
