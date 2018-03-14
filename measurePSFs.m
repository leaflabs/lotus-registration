clear all;
close all;

[~,hostname] = system('hostname')

if ~isempty(strfind(hostname, 'Justins-Mac'))
	ppath = '/Users/justin/Desktop/DDLFM';
elseif ~isempty(strfind(hostname, 'willis'))
	ppath = '/home/jkinney/Desktop/DDLFM';
else
	disp('WTF?');
    keyboard
end

timestamp = [datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%s/%s.log',ppath, timestamp);
diary(fname)
tic

ddlfm = {
    'bead/20170915',
    %'bead/20180109', 40x
    'bead/20180227',
    'bead/20180302',
    };
%%
if 0 > 1
    for i=1:numel(ddlfm)
        opath = [ppath '/' ddlfm{i} '/registration/'];
        p = [opath '*_multiply_sqrt.mat'];
        d = dir(p);
        fprintf('%d mat files found for\n%s\n\n',numel(d),p);
        getbeads(d,opath);
    end
end

if 0 < 1
config = struct;
crop_um = 20;
config.zspacing = 0.5;
config.pixel = 0.323; % um
crop_y_pixel = ceil(crop_um/config.pixel);
crop_z_pixel = ceil(crop_um/config.zspacing);
config.crop_z = crop_z_pixel;
config.crop_row = crop_y_pixel;
config.div = 4;
config.latBuf = 1;
config.axialBuf = 1;
config.threshold = 0.5;
config.figpos = [230        -895        1440         823];
config.clip = [100,100,100,100,0,0];

blacklist = [];
bead_cell_DDLFM = {};
bead_cell_LFM2 = {};
for i=1:numel(ddlfm)
    % process DDLFM
    opath = [ppath '/' ddlfm{i} '/registration/'];
    p = [opath '*_bead_xy.mat'];
    d = dir(p);
    fprintf('%d mat files found for\n%s\n\n',numel(d),p);
    for j=1:numel(d)
        [bead_psf_DDLFM, bead_psf_LFM2] = axial_psf (d(j), blacklist, config);
        bead_cell_DDLFM = [bead_cell_DDLFM {bead_psf_DDLFM}];
        bead_cell_LFM2 = [bead_cell_LFM2 {bead_psf_LFM2}];
    end
end

save('temp.mat','bead_cell_DDLFM', 'bead_cell_LFM2', 'ppath', 'timestamp');

else
    load('temp.mat');
end

plotAllPSFs(bead_cell_DDLFM, bead_cell_LFM2, ppath, timestamp);

elapsedTime = toc
diary off;

function out = plotAllPSFs (bead_cell_DDLFM, bead_cell_LFM2, outpath, timestamp)
% For reference:
%    bead_cell_LFM2 = [bead_cell_LFM2 {bead_psf_LFM2}];
%    bead_psf_DDLFM = [bead_psf_DDLFM;[z, bead_center_y_1, fwhm_y_1, 
%                                         bead_center_z, fwhm_z,
%                                         bead_center_y_2, fwhm_y_2,
%                                         bead_center_x, fwhm_x]];

bead_mat_DDLFM = cat(1,bead_cell_DDLFM{:});

z_DDLFM      = bead_mat_DDLFM(:,1);
pos_1_DDLFM  = bead_mat_DDLFM(:,2);
fwhm_1_DDLFM = bead_mat_DDLFM(:,3);
pos_3_DDLFM  = bead_mat_DDLFM(:,4);
fwhm_3_DDLFM = bead_mat_DDLFM(:,5);
pos_1alt_DDLFM  = bead_mat_DDLFM(:,6);
fwhm_1alt_DDLFM = bead_mat_DDLFM(:,7);
pos_2_DDLFM  = bead_mat_DDLFM(:,8);
fwhm_2_DDLFM = bead_mat_DDLFM(:,9);

bead_mat_LFM2 = cat(1,bead_cell_LFM2{:});

z_LFM2      = bead_mat_LFM2(:,1);
pos_1_LFM2  = bead_mat_LFM2(:,2);
fwhm_1_LFM2 = bead_mat_LFM2(:,3);
pos_3_LFM2  = bead_mat_LFM2(:,4);
fwhm_3_LFM2 = bead_mat_LFM2(:,5);
pos_1alt_LFM2  = bead_mat_LFM2(:,6);
fwhm_1alt_LFM2 = bead_mat_LFM2(:,7);
pos_2_LFM2  = bead_mat_LFM2(:,8);
fwhm_2_LFM2 = bead_mat_LFM2(:,9);

%% fwhm dim i vs pos dim i
h = figure();

subplot(3,1,1);
hold on;
plot(pos_1_DDLFM,fwhm_1_DDLFM,'k*');
plot(pos_1_LFM2,fwhm_1_LFM2,'rd');
hold off;
xlabel('dim 1 or height (um)');
ylabel('FWHM in dim 1 (um)');
title(['Full-width at Half-Max,  N = ' num2str(length(bead_mat_DDLFM)) ' beads']);

subplot(3,1,2);
hold on;
plot(pos_2_DDLFM,fwhm_2_DDLFM,'k*');
plot(pos_2_LFM2,fwhm_2_LFM2,'rd');
hold off;
xlabel('dim 2 or width (um)');
ylabel('FWHM in dim 2 (um)');

subplot(3,1,3);
hold on;
plot(pos_3_DDLFM,fwhm_3_DDLFM,'k*');
plot(pos_3_LFM2,fwhm_3_LFM2,'rd');
hold off;
xlabel('dim 3 or depth (um)');
ylabel('FWHM in dim 3 (um)');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
myStr = 'DDLFM = black asterick, LFM2 = red diamond';
text(0.4,0.9,myStr,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

filename=sprintf('%s/%s_fwhm_vs_pos.png',outpath,timestamp);
print(h,filename,'-dpng');
fprintf('\n# Output file = %s\n',filename);

%% fwhm dim i vs fwhm dim j
h = figure();

subplot(2,1,1);
daspect([1,1,1]);
hold on;
plot(fwhm_2_DDLFM,fwhm_1_DDLFM,'k*');
plot(fwhm_2_LFM2,fwhm_1_LFM2,'rd');
hold off;
xl = xlim;
yl = ylim;
ylim([0 yl(2)]);
m = min([xl(2) yl(2)]);
hold on;
plot([1 m],[1 m],'b--');
hold off;
xlabel('fwhm dim 2 or width (um)');
ylabel('fwhm dim 1 or height (um)');
title(['Full-width at Half-Max,  N = ' num2str(length(bead_mat_DDLFM)) ' beads']);

subplot(2,1,2);
daspect([1,1,1]);
hold on;
plot(fwhm_3_DDLFM,fwhm_1_DDLFM,'k*');
plot(fwhm_3_LFM2,fwhm_1_LFM2,'rd');
hold off;
xl = xlim;
yl = ylim;
ylim([0 yl(2)]);
m = min([xl(2) yl(2)]);
hold on;
plot([1 m],[1 m],'b--');
hold off;
xlabel('fwhm dim 3 or depth (um)');
ylabel('fwhm dim 1 or height (um)');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
myStr = 'DDLFM = black asterick, LFM2 = red diamond';
text(0.4,0.9,myStr,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

filename=sprintf('%s/%s_fwhm_vs_fwhm.png',outpath,timestamp);
print(h,filename,'-dpng');
fprintf('\n# Output file = %s\n',filename);
end
