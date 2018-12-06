close all;
clear all;

config = struct;
fvec = {'04-Oct-2017 17_07_02_multiply.mat'};
config.fpath = '/Users/justin/Desktop/DDLFM/fig2bc/mat_files/';
config.xypath = config.fpath;
config.outpath = '/Users/justin/Desktop/DDLFM/fig2bc/out/';
config.zspacing = 0.5;
config.xyname = 'xycoords_side2.mat';
config.crop_um = 20;

bead_cell = {};
for i=1:length(fvec)
    config.fname = fvec{i};
    bead_psf = axial_psf(config);
    bead_cell = [bead_cell {bead_psf}];
    keyboard
end

keyboard

zspacing = 4;
crop_um = 80;

fname = 'raw_side2.mat';
xyname = 'xycoords_side2.mat';
close all;
axial_psf
bead_cell = [bead_cell {bead_psf}];
fname = 'raw_side1.mat';
xyname = 'xycoords_side1.mat';
close all;
axial_psf
bead_cell = [bead_cell {bead_psf}];


%              #   #       #       #   #
glyph = ['^','*','*','.','*','+','<','>'];
colors = ['k','k','b','b','r','r','g','g'];
figure;
indfvec = [2 3 5];
ind = [indfvec 7 8];
hold on;
for i=1:length(ind)
    keyboard
    plot ( bead_cell{ind(i)}(:,4), bead_cell{ind(i)}(:,2), glyph(ind(i)), ...
        'Color', colors(ind(i)));
    %mystr = 
end
xmax = 50;
xlim([0 xmax]);
xlabel('full-width half-max in z (um)');
ylabel('full-width half-max in y (um)');
a=ylim;
%ylim([0 N]);
legend([fvec(indfvec);'';'side1';'side2'], 'Interpreter', 'none');
plot([0;a(2)],[0;a(2)],'k-.');
hold off;

i = 2;
z_mul_mean = mean(bead_cell{i}(:,4));
y_mul_mean = mean(bead_cell{i}(:,2));
n_mul_mean = length(bead_cell{i}(:,4));
i = 7;
z_raw_side1_mean = mean(bead_cell{i}(:,4));
y_raw_side1_mean = mean(bead_cell{i}(:,2));
n_raw_side1_mean = length(bead_cell{i}(:,4));
i = 8;
z_raw_side2_mean = mean(bead_cell{i}(:,4));
y_raw_side2_mean = mean(bead_cell{i}(:,2));
n_raw_side2_mean = length(bead_cell{i}(:,4));

disp(sprintf('%10s %10s %10s %10s','vol','ymean(um)','zmean(um)','N'));
disp(sprintf('%10s %10.1f %10.1f %10d','side1',y_raw_side1_mean,z_raw_side1_mean,n_raw_side1_mean));
disp(sprintf('%10s %10.1f %10.1f %10d','side2',y_raw_side2_mean,z_raw_side2_mean,n_raw_side2_mean));
disp(sprintf('%10s %10.1f %10.1f %10d','multiply',y_mul_mean,z_mul_mean,n_mul_mean));

