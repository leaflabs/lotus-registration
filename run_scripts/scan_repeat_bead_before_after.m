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

ddlfm_path = 'bead/20180328/registration/';
opath = [ppath '/' ddlfm_path];

savePath = [opath 'analysis/'];
if ~exist(savePath,'dir')
    status = mkdir(savePath);
    if status == 1
        disp(['Created folder: ' savePath]);
    else
        disp(['Error attempting to create folder:' savePath]);
        status
        exit;
    end
end

timestamp = [datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%sscan_repeat_bead_before_after_%s.log',savePath, timestamp);
diary(fname)
tic

tmp_f = 'tmp.mat';

%% Gather data from files
if ~exist([savePath tmp_f],'file')
    % % row col z
    % % y x z
    % % dim 1 2 3
    beads = {};
    out = {};
    p = [opath '*.log'];
    d = dir(p);
    fprintf('%d log files found for\n%s\n\n',numel(d),p);
    N = length(d);
    for i=1:N
        fprintf('file %d of %d\n',i,N);
        fname = [d(i).folder '/' d(i).name]
        % read log to get registration parameters for volume
        out{i} = read_log (fname);
        out{i}
    end    
    save([savePath tmp_f]);
else
    load([savePath tmp_f],'out','d');
end

%% plot
symbols = {'o','o','o','o','o','o','o','o','o','o','o','o','o'};
%symbols = {'o'};
colors = {[0 0 0]
    [1 0 0]
    [0 1 0]
    [0 0 1]
    [1 0 1]
    [0 1 1]
    [0.5 0.5 0.5]
    [0.8 0.2 0.1]};
% colors = {[0 0 0],
%     [0 0 0],
%     [0 0 0],
%     [0 0 0],
%     [0 0 0],
%     [0 0 0],
%     [0 0 0],
%     [0 0 0]};
sets = {'before_10'
    'before_9'
    'before_8'
    'before_6'
    'after_10'
    'after_9'
    'after_8'
    'after_6'};

datamat = zeros(numel(out),6);
for i=1:numel(out)
    datamat(i,:) = [out{i}.offsetF(1) out{i}.offsetF(2) out{i}.offsetF(3)...
    out{i}.rot(1) out{i}.rot(2) out{i}.rot(3)];
end
myylim = [max(datamat)
    min(datamat)];
myylim = [-100 -100 -100 -90 -90 -90;
    100 100 100 270 90 90];

plot_raw(out, savePath, colors, symbols, d, sets, myylim, timestamp);
plot_good(out, savePath, colors, symbols, d, sets, myylim, timestamp);
plot_off23(out, savePath, colors, symbols, d, sets, myylim, timestamp);

diary off;

%% Functions

%% plot offset dim 2 vs dim 3
function plot_off23 (out, ppath, colors, symbols, d, sets, myylim, timestamp)
off23 = [];
good = [1 2 4 5];
for i=1:numel(out)
    [color,ind] = getColor(d(i).name, colors);
    if isempty(find(ind==good))
        continue;
    end
    off23 = [off23;
    out{i}.offsetF(2)*0.5/0.323 out{i}.offsetF(3) ];
end

f = figure;
hold on;
x = off23(:,1);
y = off23(:,2);
plot(x,y,'o');
xlabel('offset dim 2 scaled by 0.5/0.323 [um]');
ylabel('offset dim 3 [um]');
Fit = polyfit(x,y,1); % x = x data, y = y data, 1 = order of the polynomial.
xrange = [-60 0];
yfit = polyval(Fit,xrange);
plot(xrange,yfit);
xlim(xrange);
ylim(xrange);
plot(xrange,xrange,'--');
hold off;

str = sprintf('%s/offset_dim2_vs_dim3_%s.png', ppath, timestamp);
title(str,'Interpreter','none');
ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt = sprintf('y-intercept = %2.2f',Fit(2));
text(0.7,0.5,txt,'FontSize',12,'Color',[0 0 0],'Interpreter','none');
txt = sprintf('slope = %2.2f',Fit(1));
text(0.7,0.4,txt,'FontSize',12,'Color',[0 0 0],'Interpreter','none');
print(f,str,'-dpng');

end

function out = getDataSetName (log_filename)
content = log_filename;
expr = '[^\n]*Mono';
match = regexp(content,expr,'match');
% Assume example: Recon3D_after_10_Mono
out = match{1}(9:end-5);
end

function out = getInputFileName (mat_filename)
content=mat_filename;
expr = '[^\n]*__';
match = regexp(content,expr,'match');
out = [match{1}(1:end-2) '.mat'];
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

function out = getSize (v)
out = v * v * 4;
%out = 12;
end

function [color,index] = getColor (fname, colors)
index = -1;
if ~isempty(strfind(fname, 'before_10'))
    %color = 'k';
    %color = [0 0 0];
    color = colors{1};
    index = 0;
elseif ~isempty(strfind(fname, 'before_9'))
    %color = 'r';
    %color = [1 0 0];
    color = colors{2};
    index = 1;
elseif ~isempty(strfind(fname, 'before_8'))
    %color = 'b';
    %color = [0 1 0];
    color = colors{3};
    index = 2;
elseif ~isempty(strfind(fname, 'before_6'))
    %color = 'g';
    %color = [0 0 1];
    color = colors{4};
    index = 3;
elseif ~isempty(strfind(fname, 'after_10'))
    %color = 'r';
    %color = [1 0 1];
    color = colors{5};
    index = 4;
elseif ~isempty(strfind(fname, 'after_9'))
    %color = 'r';
    %color = [0 1 1];
    color = colors{6};
    index = 5;
elseif ~isempty(strfind(fname, 'after_8'))
    %color = 'r';
    %color = [0.5 0.5 0.5];
    color = colors{7};
    index = 6;
elseif ~isempty(strfind(fname, 'after_6'))
    %color = 'r';
    %color = [1 1 0];
    color = colors{8};
    index = 7;
else
    color = [0 0 0];
end
end

function title = plotOffset (out, symbols, d, colors, index, indices)
hold on;
for i=1:numel(out)
    % set color
    %color = getColor(d(i).name, colors{i});
    [color,ind] = getColor(d(i).name, colors);
    if isempty(find(ind==indices))
        continue;
    end
    % set symbol
    %symbol = symbols{i};
    symbol = 'o';
    % set symbol size
    %s = getSize(out{i}.final_MI_frac);
    s = 5;
    % plot
    if index<3
        plot( ind, out{i}.offsetF(index)*0.5/0.323, 'Marker', symbol, 'Color', color, 'MarkerSize',s);
        %plot( ind, out{i}.offsetI(index), 'Marker', '.', 'Color', [1 0 0], 'MarkerSize',8);
    else
        plot( ind, out{i}.offsetF(index), 'Marker', symbol, 'Color', color, 'MarkerSize',s);
    end
end
hold off;
title = 'before = BLACK, after = RED';

end

function plotRot (out, symbols, d, colors, index, indices)
hold on;
for i=1:numel(out)
    % set color
    %color = getColor(d(i).name, colors{i});
    [color,ind] = getColor(d(i).name, colors);
    if isempty(find(ind==indices))
        continue;
    end
    % set symbol
    symbol = 'o';
    %symbol = symbols{j};
    % set symbol size
    %s = getSize(out{i}.final_MI_frac);
    s = 5;
    % plot
    plot( ind, out{i}.rot(index)*180/pi, 'Marker', symbol, 'Color', color, 'MarkerSize',s);
    plot( ind, out{i}.angle(index)*180/pi, 'Marker', '.', 'Color', [0.2 0.4 0.2], 'MarkerSize',8);
end
hold off;
end




%%
function plot_raw (out, ppath, colors, symbols, d, sets, myylim, timestamp)

%a = numel(out);
f = figure;

good = [0:7];

subplot(1,3,1);
plotOffset (out, symbols, d, colors, 1, good);
ylabel('offset in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
mytitle = plotOffset (out, symbols, d, colors, 2, good);
ylabel('offset in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end
%title(mytitle);
str = [ppath '/distribution_offset.png'];
title(str,'Interpreter','none');

subplot(1,3,3);
plotOffset (out, symbols, d, colors, 3, good);
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
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,1));

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,2));

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,3));

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% for i=1:numel(sets)
% text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
% end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
for i=1:numel(sets)
    text(0.165+(i-1)*0.0175,0.01,sets{i},'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
end

str = sprintf('%s/distribution_offset_%s.png', ppath, timestamp);
print(f,str,'-dpng');

% rotation
h = figure;

subplot(1,3,1);
plotRot (out, symbols, d, colors, 1, good);
ylabel('rotation around dim 1 (degrees)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
plotRot (out, symbols, d, colors, 2, good);
ylabel('rotation around dim 2 (degrees)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end
%title(mytitle);
str = [ppath '/distribution_rotation.png'];
title(str,'Interpreter','none');

subplot(1,3,3);
plotRot (out, symbols, d, colors, 3, good);
ylabel('rotation around dim 3 (degrees)');
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
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,4));

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,5));

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,6));

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% for i=1:numel(sets)
% text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
% end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
for i=1:numel(sets)
    text(0.165+(i-1)*0.0175,0.01,sets{i},'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
end

str = sprintf('%s/distribution_rotation_%s.png', ppath, timestamp);
print(h,str,'-dpng');

end

%%
function plot_good (out, ppath, colors, symbols, d, sets, myylim, timestamp)

%a = numel(out);
f = figure;

good = [1 2 4 5];

subplot(1,3,1);
plotOffset (out, symbols, d, colors, 1, good);
ylabel('offset in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
mytitle = plotOffset (out, symbols, d, colors, 2, good);
ylabel('offset in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end
%title(mytitle);
str = [ppath '/distribution_offset.png'];
title(str,'Interpreter','none');

subplot(1,3,3);
plotOffset (out, symbols, d, colors, 3, good);
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
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,1));

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,2));

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,3));

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% for i=1:numel(sets)
% text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
% end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
for i=1:numel(sets)
    text(0.165+(i-1)*0.0175,0.01,sets{i},'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
end

str = sprintf('%s/distribution_offset_good_%s.png', ppath, timestamp);
print(f,str,'-dpng');

% rotation
h = figure;

subplot(1,3,1);
plotRot (out, symbols, d, colors, 1, good);
ylabel('rotation around dim 1 (degrees)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
plotRot (out, symbols, d, colors, 2, good);
ylabel('rotation around dim 2 (degrees)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end
%title(mytitle);
str = [ppath '/distribution_rotation.png'];
title(str,'Interpreter','none');

subplot(1,3,3);
plotRot (out, symbols, d, colors, 3, good);
ylabel('rotation around dim 3 (degrees)');
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
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,4));

subplot(1,3,2);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,5));

subplot(1,3,3);
y = ylim;
ndel = y(2)-y(1);
m = (del-ndel)/2;
ylim([y(1)-m y(2)+m]);
x = xlim;
xlim([x(1)-1 x(2)]);
ylim(myylim(:,6));

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% for i=1:numel(sets)
% text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
% end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
for i=1:numel(sets)
    text(0.165+(i-1)*0.0175,0.01,sets{i},'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
end

str = sprintf('%s/distribution_rotation_good_%s.png', ppath, timestamp);
print(h,str,'-dpng');

end



%%
% function calc_stats (out, ppath, colors, d, sets)
% 
% 
% 
% % for each data set
% 
% good = [0 1 2 4 5];
% 
% for i=1:numel(out)
%     [color,ind] = getColor(d(i).name, colors);
%     if isempty(find(ind==good))
%         continue;
%     end
%     % set symbol
%     vec{ind}. out{i}.offsetF(index), 'Marker', symbol, 'Color', color, 'MarkerSize',s);
%     %plot( ind, out{i}.offsetI(index), 'Marker', '.', 'Color', [1 0 0], 'MarkerSize',8);
% end
% 
% 
% subplot(1,3,1);
% plotOffset (out, symbols, d, colors, 1, good);
% ylabel('offset in dim 1 (um)');
% set(gca,'xtick',[]);
% y = ylim;
% del = y(2)-y(1);
% 
% subplot(1,3,2);
% mytitle = plotOffset (out, symbols, d, colors, 2, good);
% ylabel('offset in dim 2 (um)');
% set(gca,'xtick',[]);
% y = ylim;
% ndel = y(2)-y(1);
% if ndel > del
%     del = ndel;
% end
% %title(mytitle);
% str = [ppath '/distribution_offset.png'];
% title(str,'Interpreter','none');
% 
% subplot(1,3,3);
% plotOffset (out, symbols, d, colors, 3, good);
% ylabel('offset in dim 3 (um)');
% set(gca,'xtick',[]);
% y = ylim;
% ndel = y(2)-y(1);
% if ndel > del
%     del = ndel;
% end
% 
% subplot(1,3,1);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% x = xlim;
% xlim([x(1)-1 x(2)]);
% 
% subplot(1,3,2);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% x = xlim;
% xlim([x(1)-1 x(2)]);
% 
% subplot(1,3,3);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% x = xlim;
% xlim([x(1)-1 x(2)]);
% 
% % ax1 = axes('Position',[0 0 1 1],'Visible','off');
% % axes(ax1);
% % for i=1:numel(sets)
% % text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
% % end
% 
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% for i=1:numel(sets)
%     text(0.165+(i-1)*0.0175,0.01,sets{i},'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
% end
% 
% str = [ppath '/distribution_offset_good.png'];
% print(f,str,'-dpng');
% 
% % rotation
% h = figure;
% 
% subplot(1,3,1);
% plotRot (out, symbols, d, colors, 1, good);
% ylabel('rotation around dim 1 (degrees)');
% set(gca,'xtick',[]);
% y = ylim;
% del = y(2)-y(1);
% 
% subplot(1,3,2);
% plotRot (out, symbols, d, colors, 2, good);
% ylabel('rotation around dim 2 (degrees)');
% set(gca,'xtick',[]);
% y = ylim;
% ndel = y(2)-y(1);
% if ndel > del
%     del = ndel;
% end
% %title(mytitle);
% str = [ppath '/distribution_rotation.png'];
% title(str,'Interpreter','none');
% 
% subplot(1,3,3);
% plotRot (out, symbols, d, colors, 3, good);
% ylabel('rotation around dim 3 (degrees)');
% set(gca,'xtick',[]);
% y = ylim;
% ndel = y(2)-y(1);
% if ndel > del
%     del = ndel;
% end
% 
% subplot(1,3,1);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% x = xlim;
% xlim([x(1)-1 x(2)]);
% 
% subplot(1,3,2);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% x = xlim;
% xlim([x(1)-1 x(2)]);
% 
% subplot(1,3,3);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% x = xlim;
% xlim([x(1)-1 x(2)]);
% 
% % ax1 = axes('Position',[0 0 1 1],'Visible','off');
% % axes(ax1);
% % for i=1:numel(sets)
% % text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
% % end
% 
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% for i=1:numel(sets)
%     text(0.165+(i-1)*0.0175,0.01,sets{i},'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
% end
% 
% str = [ppath '/distribution_rotation_good.png'];
% print(h,str,'-dpng');
% 
% end


%%
function plot_beads (out, beads, ppath, colors, voxel)

a = numel(beads);
h = figure;

subplot(1,3,1);
hold on;
for i=1:a
    b = size(beads{i});
    plot( 0, beads{i}(:,1)*voxel(1), [colors{1} 'o']);
end
hold off;
ylabel('bead position in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
hold on;
for i=1:a
    b = size(beads{i});
    plot( 0, beads{i}(:,2)*voxel(2), [colors{1} 'o'] );
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
    b = size(beads{i});
    plot( 0, beads{i}(:,3)*voxel(3), [colors{1} 'o'] );
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

% %%
% h = figure;
% 
% subplot(1,3,1);
% hold on;
% for i=1:a
%     b = size(spread{i});
%     plot( zeros(b(1),1), spread{i}(:,1)./fwhm{i}(:,1), [colors{i} 'o'] );
% end
% hold off;
% ylabel('spread of bead position in dim 1 (-)');
% set(gca,'xtick',[]);
% y = ylim;
% del = y(2)-y(1);
% 
% subplot(1,3,2);
% hold on;
% for i=1:a
%     b = size(spread{i});
%     plot( zeros(b(1),1), spread{i}(:,2)./fwhm{i}(:,2), [colors{i} 'o'] );
% end
% hold off;
% ylabel('spread of bead position in dim 2 (-)');
% set(gca,'xtick',[]);
% y = ylim;
% ndel = y(2)-y(1);
% if ndel > del
%     del = ndel;
% end
% 
% subplot(1,3,3);
% hold on;
% for i=1:a
%     b = size(spread{i});
%     plot( zeros(b(1),1), spread{i}(:,3)./fwhm{i}(:,3), [colors{i} 'o'] );
% end
% hold off;
% ylabel('spread of bead position in dim 3 (-)');
% set(gca,'xtick',[]);
% y = ylim;
% ndel = y(2)-y(1);
% if ndel > del
%     del = ndel;
% end
% 
% subplot(1,3,1);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% 
% subplot(1,3,2);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% 
% subplot(1,3,3);
% y = ylim;
% ndel = y(2)-y(1);
% m = (del-ndel)/2;
% ylim([y(1)-m y(2)+m]);
% 
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% str = sprintf('Spread of bead position normalized by fwhm (N = %d bead samples)',a);
% text(0.15,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
% str = sprintf('Different colors are different samples.');
% text(0.15,0.06,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
% 
% str = [ppath '/normalized_spread_of_bead_position.png'];
% print(h,str,'-dpng');



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
eval(['out.final_MI = ' get_match_MI('Final MI',content) ';']);
eval(['out.final_MI_frac = ' get_match_MI_frac('MI frac',content) ';']);
end

% function out = read_log_MI (f)
% out = [];
% content = fileread(f);
% keyboard
% eval(['out.trans = ' get_match('trans',content) ';']);
% eval(['out.angle = ' get_match('angle',content) ';']);
% eval(['out.centroid = ' get_match('centroid',content) ';']);
% eval(['out.offsetI = ' get_match('offset',content) ';']);
% out.offsetF = out.trans - out.centroid;
% eval(['out.rot = ' get_match('rot',content) ';']);
% end

function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{1}{1};
end

function out = get_match_MI (str, content)
expr = [str ' = [^\n]* final offset'];
match = regexp(content,expr,'match');
C = strsplit(match{end},' ');
out = single(C{4});
end

function out = get_match_MI_frac (str, content)
expr = [str ' = [^\n]*,'];
match = regexp(content,expr,'match');
C = strsplit(match{end},' ');
out = single(C{4}(1:end-1));
end