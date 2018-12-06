clear all;
close all;

%fname = '106_Recon3D_20170322-CROPPED_3x3_Nnum_15__27-Mar-2017_14:01:33.mat';
fname = 'raw_side2.mat';
%fname = 'raw_side2.mat';
%fpath = '/Users/justin/Desktop/lotus/lotus-data/';
fpath = '/home/jkinney/Desktop/LFM volume registration/from_nikita/';

%microlens_array_pitch = 500.0; % um
%mag = 50.0;
%Nnum = 15;
%zspacing = 10;
%maxIter = 8;
%zslices = 10;
%zrange = [-300, 300];

load([fpath fname]);
xy = squeeze(sum(XguessSAVE1,3));
figure(1);
imagesc(xy,[0 2^8]);
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
save('xycoords_side2.mat');

% Note that I accidentally closed the figure which terminated this program
% before saving the mat file. Fortunately, I was able to use the stdout
% data to create user_click.m, in which I save to mat file.