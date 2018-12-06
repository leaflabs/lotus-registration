function bead_psf = axial_psf (config)

timestamp = datestr(datetime('now'));
fdiary = sprintf('%s%s.log',config.outpath,timestamp);
%fdiary = 'axial_psf.log'
diary(fdiary)
disp('%%');
disp(['%% ' datestr(datetime)]);
disp('%%');

%fname = '106_Recon3D_20170322-CROPPED_3x3_Nnum_15__27-Mar-2017_14:01:33.mat';
%fpath = '/Users/justin/Desktop/lotus/lotus-data/';


% a = who;
% keyboard
% for i=1:length(a)
%     mystr = sprintf('config.%s = %s;',a{i},a{i});
%     disp(mystr);
%     eval(mystr)
% end
%config.fname = fname;
%config.fpath = fpath;
%config.fname = 'raw_side1.mat';
%config.fpath = '/home/jkinney/Desktop/LFM volume registration/from_nikita/';
%config.xyname = 'xycoords.mat';
%config.xyname = xyname;
%config.crop_z = crop_z;
%config.xypath = '/home/jkinney/Desktop/LFM volume registration/psf/';
%config.outpath = outpath;
%config.PSFpath = '/home/jkinney/Desktop/lotus/lotus-recon/PSFmatrix/';
%config.PSFname = 'PSFmatrix_Z_from_-300_to_300_zspacing_10_Nnum_15_OSR_3__24-Mar-2017_16:11:44.mat';

%config.microlens_array_pitch = 500.0; % um
config.pixel = 0.323; % um
%config.mag = 50.0;
%config.Nnum = 15;
%config.zspacing = 4; % um
%config.zspacing = zspacing; % um
%config.maxIter = 8;
config.zslices = 60;
config.zrange = [-300, 300];
% JPK
%config.crop_um = crop_um; % um
%config.crop_um = 60; % um
config.crop_y_pixel = ceil(config.crop_um/config.pixel);
config.crop_z_pixel = ceil(config.crop_um/config.zspacing);
config.crop_z = config.crop_z_pixel;
config.crop_row = config.crop_y_pixel;
config.div = 4;
%JPK
%config.crop_row = 101; % KEEP ODD, size in pixels of cropped image 
%config.crop_z = 11; % KEEP ODD, size in pixels of cropped image
config.latBuf = 1;
config.axialBuf = 1;
config.threshold = 0.5;
config.figpos = [230        -895        1440         823];

%blacklist = [13 15 16 18 25 31 34 35 37 44 45 53 66 67];

blacklist = [];

tmp_f = [config.outpath 'temp.mat'];

if ~exist(tmp_f,'file')
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
        if ~isempty(find(blacklist==i))
            msg = sprintf('Skipping index %d since blacklisted.',i);
            disp(msg);
            continue;
        end
        % find z value of bead
        row = xycoords(i,2);
        col = xycoords(i,1);
        z = [z get_z( row, col, XguessSAVE1 )];
        if z(end) < 0
            disp('Error! No z value found.');
            continue
        end
        disp(sprintf('\n\n'));
        msg = sprintf('bead %d found at row = %d, col = %d, z = %d',length(index),row,col,z(end));
        disp(msg);
        %
        % slice volume
        csec = squeeze(XguessSAVE1(:,col,:));
        %
        latBuf = 1;
        axialBuf = 1;
        % = crop_index = [minrow maxrow mincol maxcol]
        crop_index = cropIndex(config, size(csec), row, z(end));
        if crop_index<0
            continue;
        end
        csec_cropped = cropImag(csec, crop_index);
        % plot image
        h = figure;
        %set(gcf,'Position',config.figpos);
        [dim, ztick_values, ztick_labels] = plot_slice_cropped ( config, h, ...
            csec_cropped, size(csec), crop_index);
%          [pixel,allLabels,dim,nz] = plot_slice_cropped ( config, h, ...
%             csec_cropped, size(csec), crop_index);
        lateralVals = plot_lat_psf_cropped(config, h, csec_cropped);
        axialVals = plot_axial_psf_cropped(config, h, csec_cropped);
        % display surface
        %subplot(2,2,4);
        %surf(csec);
        %title(msg);
        %
        % calculate PSFs
        latC = latPSF (config, h, lateralVals);
        %axialC = axialPSF (config, h, axialVals,allLabels,dim,nz,size(csec_cropped),z(end));
        axialC = axialPSF (config, h, axialVals, crop_index,...
            dim,size(csec_cropped),z(end),0);
        %
        if latC(2)>0 & axialC(2)>0
            bead_psf = [bead_psf;[latC(1) latC(2) axialC(1) axialC(2)]];
            index = [index i];
            % save figure
            %str=sprintf('%s%s_bead%03d.pdf',fpath,fname(1:end-3),i);
            %print(gcf,str,'-dpdf');
            %disp(sprintf('# Output file = %s',str));
        end
        figure(h);
        text(50.0,0.5,sprintf('bead %d',i),'FontSize',12);
        %
        str=sprintf('%s%s_bead%03d.png',config.outpath,config.fname(1:end-4),i);
        f=getframe(gcf);
        [X, map] = frame2im(f);
        imwrite(X, str);
        keyboard
    end
    save(tmp_f);
else
    load(tmp_f);
end

keyboard

h = figure();
plotAllPSFs(h, bead_psf, config.fpath);
str=sprintf('%s%s_all_beads.png',config.outpath,config.fname(1:end-4));
f=getframe(gcf);
[X, map] = frame2im(f);
imwrite(X, str);
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
msg = sprintf('Orig slice is %d row by %d z',a(1),a(2));
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


function axialVals = plot_axial_psf_cropped (config, handle, csec)
subplot(2,2,3);
a = max(csec);
axialVals = single(a)/max(single(a(config.axialBuf:length(a)-config.axialBuf+1)));
plot(axialVals,'*-');
set(gca, 'YDir', 'reverse');
%thisdim = get(gca,'Position');
%set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
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

function [dim,z,Lzc] = plot_slice_cropped ( config, handle, ...
    csec_cropped, a, crop_index )
figure(handle);
h = subplot(2,2,1);
im = imagesc(csec_cropped,[0 2^8]);
%colorbar();
% h = subplot(2,2,3);
% im = imagesc(csec_cropped,[0 2^8]);
%colorbar();
%
% mark bead
% BROKEN
% hold on;
% plot(z,row,'ow');
% hold off;
%
% transform figure so that plot axes have same scale
% to do so,
% the target crop size is config.crop_um, e.g. 100 um
% given crop_index calculate actual size in um
% get size of plot (in whatever units)
% scale plot so 1 um in x is same size on screen as 1 um in y
%
b = size(csec_cropped);
%pixel = config.microlens_array_pitch / single(config.Nnum) / config.mag;
pixel = config.pixel;
vertical_microns = single(b(1)) * pixel;
depth_microns = (single(b(2))-1.0) * single(config.zspacing);
% JPK
%vertical_vox = a(1);
%depth_vox = a(2);
%
disp(['image size = ' num2str(vertical_microns) ' um tall x ' num2str(depth_microns) ' um wide' ]);
% pos = [left bottom width height]
% here, width is depth, and height is vertical
pos = get(h,'Position'); %  gives x left, y bottom, width, height
newpos = pos;
% scale width
newpos(3) = pos(3) * vertical_microns / depth_microns;
%newpos(3) = pos(3) * vertical_vox / depth_vox;
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
% crop_ind = [minrow maxrow mincol maxcol]
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
%{num2str(L(1)),num2str(L(2)),num2str(L(3)),num2str(L(4)),num2str(L(5)) };

% nz = config.zslices;
% dz = config.zspacing;
% allLabels = [-nz*dz:dz:nz*dz];
%allLabels = [config.zrange(1):dz:config.zrange(2)];
%step = 3;
%tmp = [step:step:nz];
%ax.XTick = [-fliplr(tmp) 0 tmp] + (nz+1)*ones(1,2*length(tmp)+1);
% TODO FIX HARDCODING
% ax.XTick = [1:11];
% labelindex = [z-5:z+5];
% if isempty(find(labelindex>0))
%     disp('WTF?')
%     keyboard
% end
% if z==65
%     keyboard
% end
% labelindex = labelindex(find(labelindex>0));
% ax.XTickLabel = {allLabels(labelindex)};
xlabel(' z, i.e. depth (um)');
%ax.YTickLabel = {round(ax.YTick*pixel)};
ylabel(' y, i.e. height (um)');
dim = get(gca, 'Position');
% add title
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
