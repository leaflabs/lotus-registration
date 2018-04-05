clear all;
close all;

index_pattern = 'Mono';
%ppath = '/home/jkinney/Desktop/DDLFM/bead/20180328/registration';
ppath = '/Users/justin/Desktop/DDLFM/bead/20180328/registration';
%ppath = '/home/jkinney/Desktop/DDLFM/bead/20180319/registration';
pathNpattern = sprintf('%s/*before*.log',ppath);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end

% /home/jkinney/Desktop/DDLFM/bead/20180319/registration/Recon3D_1_Mono_N15__20180321_173750_multiply_sqrt.mat
% DDLFM: bead 1 found at row = 258, col = 342, z = 269
% DDLFM: bead 2 found at row = 718, col = 244, z = 446
% DDLFM: bead 3 found at row = 1012, col = 443, z = 319
% /home/jkinney/Desktop/DDLFM/bead/20180319/registration/Recon3D_2_Mono_N15__20180321_174905_multiply_sqrt.mat
% DDLFM: bead 1 found at row = 956, col = 357, z = 386
% DDLFM: bead 2 found at row = 490, col = 500, z = 259
% /home/jkinney/Desktop/DDLFM/bead/20180319/registration/Recon3D_3_Mono_N15__20180321_180959_multiply_sqrt.mat
% DDLFM: bead 1 found at row = 751, col = 273, z = 527
% DDLFM: bead 2 found at row = 718, col = 382, z = 376
% DDLFM: bead 3 found at row = 494, col = 441, z = 364
% /home/jkinney/Desktop/DDLFM/bead/20180319/registration/Recon3D_4_Mono_N15__20180321_182118_multiply_sqrt.mat
% DDLFM: bead 1 found at row = 605, col = 213, z = 520
% DDLFM: bead 2 found at row = 1005, col = 97, z = 658
% DDLFM: bead 3 found at row = 735, col = 365, z = 367

% % row col z
% % y x z
% % dim 1 2 3
% beads = {
% [258,  342, 269,
% 718,  244, 446,
% 1012, 443, 319]
%  
% [956, 357, 386,
% 490, 500, 259]
%  
% [  751,  273, 527,
%   718,  382, 376,
%   494,  441, 364]
%  
% [ 605,  213, 520,
%  1005,  97, 658,
%  735, 365, 367]
% };
% 
% fwhm = {
% [ 8 9 12 
% 8 10 12 
% 10 10 12 ]
% [ 7 6 10
% 9 10 11 ]
% [ 8 10 11
% 8 8 9 
% 10 10 10 ]
% [ 8 10 11 
% 8 10 11 
% 9 10 11]
% };

if 0<1
    o = {[],[],[],[],[],[],[],[],[],[]};
    r = {[],[],[],[],[],[],[],[],[],[]};
    offset = {};
    bl = {[],[],[],[],[],[],[],[],[],[]};
    %angle = {}
    n = numel(logs);
    for i=1:n
        f = [logs(i).folder '/' logs(i).name];
        fprintf('%d of %d, %s\n',i,n,f);
        out = read_log (f);
        j = getIndex(logs(i).name, index_pattern);
        if j<1 & j>4
            disp('WTF?');
            keyboard
        end
        o{j} = [o{j};out.offsetF];
        r{j} = [r{j};out.rot];
        offset{j} = out.offsetI;
        % for each bead in volume
        a = size(beads{j});
        for k = 1:a(1)
            % transform bead location
            b = transformBead (beads{j}(k,:),out.centroid, out.offsetF, out.rot);
            bl{j} = [bl{j};b];
        end
    end
    save('tmp.mat');
else
    load('tmp.mat');
end

voxel_size = [0.323 0.323 4/8]; % um
spread = calcSpread(bl, beads, voxel_size);

colors = {'k','r','b','g'};
plot_raw(o, r, offset, out.angle, ppath, colors);
plot_beads(bl, ppath, colors, voxel_size );
plot_spread(spread, ppath, colors, fwhm);

%%

function spread = calcSpread (bl, beads, voxel_size)
%spread = zeros(11,3);
spread = {[],[],[],[]};

% for each sample
for i = 1:numel(bl)
    a = size(beads{i});
    n = a(1); % number of beads
    b = size(bl{i});
    % for each bead
    for j = 1:n
        repeats = bl{i}(j:n:b(1),:);
        pixel_spread = max(repeats) - min(repeats);
        um_spread = pixel_spread .* voxel_size;
        spread{i} = [spread{i};um_spread];
    end
end
end


function out = transformBead (bead, centroid, offset, rot)
a = bead - centroid;
b = rotate (a, rot); % CHECK
out = b + centroid + offset;
end

function out = rotation_matrix (angle)
a = [1 0 0;...
    0 cos(angle(1)) -sin(angle(1));...
    0 sin(angle(1)) cos(angle(1))];
b = [cos(angle(2)) 0 -sin(angle(2));...
    0 1 0;...
    sin(angle(2)) 0 cos(angle(2))];
c = [cos(angle(3)) -sin(angle(3)) 0;...
    sin(angle(3)) cos(angle(3)) 0;...
    0 0 1];
out = a*b*c;
end

function [out] = rotate (pos, angle)
rot = single( rotation_matrix (angle) );
out = pos*rot;
end

function out = getIndex (fname, pattern)
% fname = Recon3D_1_Mono_N15__20180321_173750.log
i = strfind(fname,pattern);
if strcmp(fname(i-3),'_')
    out = str2num(fname(i-2));
else
    out = str2num(fname(i-3:i-2));
end
end

%%
function plot_raw (o,r, offset, angle, ppath, colors)

a = numel(o);
h = figure;

subplot(1,3,1);
hold on;
for i=1:a
    b = size(o{i});
    plot( zeros(b(1),1), o{i}(:,1), [colors{i} 'o'] );
    plot( 0, offset{i}(1),[colors{i} '*']);
end
hold off;
ylabel('offset in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
hold on;
for i=1:a
    b = size(o{i});
    plot( zeros(b(1),1), o{i}(:,2), [colors{i} 'o'] );
    plot( 0, offset{i}(2),[colors{i} '*']);
end
hold off;
ylabel('offset in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
hold on;
for i=1:a
    b = size(o{i});
    plot( zeros(b(1),1), o{i}(:,3), [colors{i} 'o'] );
    plot( 0, offset{i}(3),[colors{i} '*']);
end
hold off;
ylabel('offset in dim 3 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,1);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Repeatability of registration parameters (N = %d bead samples)',a);
text(0.15,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
str = sprintf('Circles are final values. Asterisks are initial guesses.\nDifferent colors are different samples. Each sample was registered multiple times.');
text(0.15,0.06,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_offset.png'];
print(h,str,'-dpng');

% rotation
h = figure;

subplot(1,3,1);
% plot( zeros(a(1),1), r(:,1), 'ko' );
% hold on;
% plot( 0, out.angle(1),'ro');
% hold off;
hold on;
for i=1:a
    b = size(o{i});
    plot( zeros(b(1),1), r{i}(:,1), [colors{i} 'o'] );
    %plot( 0, angle{i}(1),[colors{i} '*']);
end
plot( 0, angle(1),'c*');
hold off;
ylabel('rotation around dim 1 (radians)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);


subplot(1,3,2);
% plot( zeros(a(1),1), r(:,2), 'ko' );
% hold on;
% plot( 0, out.angle(2),'ro');
% hold off;
hold on;
for i=1:a
    b = size(o{i});
    plot( zeros(b(1),1), r{i}(:,2), [colors{i} 'o'] );
    %plot( 0, angle{i}(1),[colors{i} '*']);
end
plot( 0, angle(2),'c*');
hold off;
ylabel('rotation around dim 2 (radians)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
% plot( zeros(a(1),1), r(:,3), 'ko' );
% hold on;
% plot( 0, out.angle(3),'ro');
% hold off;
hold on;
for i=1:a
    b = size(o{i});
    plot( zeros(b(1),1), r{i}(:,3), [colors{i} 'o'] );
    %plot( 0, angle{i}(1),[colors{i} '*']);
end
plot( 0, angle(3),'c*');
hold off;
ylabel('rotation around dim 3 (radians)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,1);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Repeatability of registration parameters (N = %d bead samples)',a);
text(0.15,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
str = sprintf('Circles are final values. Asterisks are initial guesses.\nDifferent colors are different samples. Each sample was registered multiple times.');
text(0.15,0.06,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_rotation.png'];
print(h,str,'-dpng');

end

%%
function plot_beads (bl, ppath, colors, voxel)

a = numel(bl);
h = figure;

subplot(1,3,1);
hold on;
for i=1:a
    b = size(bl{i});
    plot( zeros(b(1),1), bl{i}(:,1)*voxel(1), [colors{i} 'o'] );
end
hold off;
ylabel('bead position in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
hold on;
for i=1:a
    b = size(bl{i});
    plot( zeros(b(1),1), bl{i}(:,2)*voxel(2), [colors{i} 'o'] );
end
hold off;
ylabel('bead position in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
hold on;
for i=1:a
    b = size(bl{i});
    plot( zeros(b(1),1), bl{i}(:,3)*voxel(3), [colors{i} 'o'] );
end
hold off;
ylabel('bead position in dim 3 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,1);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Repeatability of bead position (N = %d bead samples)',a);
text(0.15,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
str = sprintf('Different colors are different samples.');
text(0.15,0.06,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_bead_positions.png'];
print(h,str,'-dpng');

end

%%
function plot_spread (spread, ppath, colors, fwhm)

a = numel(spread);
h = figure;

subplot(1,3,1);
hold on;
for i=1:a
    b = size(spread{i});
    plot( zeros(b(1),1), spread{i}(:,1), [colors{i} 'o'] );
end
hold off;
ylabel('spread of bead position in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
hold on;
for i=1:a
    b = size(spread{i});
    plot( zeros(b(1),1), spread{i}(:,2), [colors{i} 'o'] );
end
hold off;
ylabel('spread of bead position in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
hold on;
for i=1:a
    b = size(spread{i});
    plot( zeros(b(1),1), spread{i}(:,3), [colors{i} 'o'] );
end
hold off;
ylabel('spread of bead position in dim 3 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,1);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Spread of bead position (N = %d bead samples)',a);
text(0.15,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
str = sprintf('Different colors are different samples.');
text(0.15,0.06,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/spread_of_bead_position.png'];
print(h,str,'-dpng');

%%
h = figure;

subplot(1,3,1);
hold on;
for i=1:a
    b = size(spread{i});
    plot( zeros(b(1),1), spread{i}(:,1)./fwhm{i}(:,1), [colors{i} 'o'] );
end
hold off;
ylabel('spread of bead position in dim 1 (-)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
hold on;
for i=1:a
    b = size(spread{i});
    plot( zeros(b(1),1), spread{i}(:,2)./fwhm{i}(:,2), [colors{i} 'o'] );
end
hold off;
ylabel('spread of bead position in dim 2 (-)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
hold on;
for i=1:a
    b = size(spread{i});
    plot( zeros(b(1),1), spread{i}(:,3)./fwhm{i}(:,3), [colors{i} 'o'] );
end
hold off;
ylabel('spread of bead position in dim 3 (-)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,1);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Spread of bead position normalized by fwhm (N = %d bead samples)',a);
text(0.15,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
str = sprintf('Different colors are different samples.');
text(0.15,0.06,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/normalized_spread_of_bead_position.png'];
print(h,str,'-dpng');



end



%% read log file
function out = read_log (f)
out = [];

content = fileread(f);

eval(['out.trans = ' get_match('trans',content) ';']);
eval(['out.angle = ' get_match('angle',content) ';']);
eval(['out.centroid = ' get_match('centroid',content) ';']);
eval(['out.offsetI = ' get_match('offset',content) ';']);
out.offsetF = out.trans - out.centroid;
eval(['out.rot = ' get_match('rot',content) ';']);
end

function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{1}{1};
end
