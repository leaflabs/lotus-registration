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

keyboard

save('temp.mat','bead_cell_DDLFM', 'bead_cell_LFM2', 'ppath', 'timestamp');

plotAllPSFs(bead_cell_DDLFM, bead_cell_LFM2, [ppath '/psf'], timestamp);

elapsedTime = toc
diary off;

function out = plotAllPSFs (bead_cell_DDLFM, bead_cell_LFM2, outpath, timestamp)
h = figure();
subplot(3,1,2);
plot(bead_center(:,3),bead_center(:,4),'o');
xlabel('z, i.e. depth (um)');
ylabel('full-width half-max in z (um)');
text(0, 70,['N = ' num2str(length(bead_center)) ' beads'],'FontSize',10);
subplot(3,1,3);
plot(bead_center(:,3),bead_center(:,2),'o');
xlabel('z, i.e. depth (um)');
ylabel('full-width half-max in y (um)');
subplot(3,1,1);
plot(bead_center(:,4),bead_center(:,2),'o');
xlabel('full-width half-max in z (um)');
ylabel('full-width half-max in y (um)');
filename=sprintf('%s_all_beads.png',rname);
print(h,filename,'-dpng');
disp(sprintf('# Output file = %s',str));
end
