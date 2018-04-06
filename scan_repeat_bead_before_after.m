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

voxel_size = [0.323 0.323 4/8]; % um
ddlfm_path = 'bead/20180328/registration/';
opath = [ppath '/' ddlfm_path];

timestamp = [datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%s/scan_repeat_bead_before_after_%s.log',opath, timestamp);
diary(fname)
tic

%% Gather data from files
if 0>1
    % % row col z
    % % y x z
    % % dim 1 2 3
    beads = {};
    out = {};
    p = [opath '*_multiply_sqrt.mat'];
    d = dir(p);
    fprintf('%d mat files found for\n%s\n\n',numel(d),p);
    for i=1:length(d)
        fname = [d(i).folder '/' d(i).name]
        % read log to get registration parameters for volume
        lf = [fname(1:end-18) '.log'];
        out{i} = read_log (lf);
        % load volume
        load(fname);
        if exist('Xvolume','var')
            vol = Xvolume;
        elseif exist('XguessSAVE1','var')
            vol = XguessSAVE1;
        else
            disp('WTF?! Unknown data name.');
            keyboard
        end
        % find bead locations in volume
        xyz = getXYZbead(vol, opath, [d(i).name(1:end-4) '_' timestamp]);
        beads{i} = xyz;
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
        %
    end
    
    %% transform bead location
    beads_trans = {};
    % for each volume
    for i=1:numel(beads)
        beads_trans{i} = [];
        a = size(beads{i});
        % for each bead in volume
        for k = 1:a(1)
            % transform bead location
            b = transformBead (beads{i}(k,:),out{i}.centroid, out{i}.offsetF, out{i}.rot);
            beads_trans{i} = [beads_trans{i};b];
        end
    end
    
    %% Gather more data from files
    for i=1:length(d)
        fname = [d(i).folder '/' d(i).name]
        % read log to get registration parameters for volume
        lf = [fname(1:end-18) '.log'];
        out{i} = read_log (lf);
    end
    
    save([opath '/tmp.mat']);
else
    load([opath '/tmp.mat']);
end

%% plot
%symbols = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
symbols = {'o'};
colors = {};
for i=1:numel(d)
    colors{i}=rand(1,3);
end
plot_raw(out, opath, colors, symbols, d);
keyboard
plot_beads(out, beads_trans, opath, colors, voxel_size );
fwhm = [];
spread = calcSpread(bl, beads, voxel_size);
plot_spread(spread, ppath, colors, fwhm);




%% Functions

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

function out = getSize (v)
out = v * v * 4;
%out = 12;
end

function color = getColor (fname, C)
if 0<1
    if ~isempty(strfind(fname, 'before'))
        color = 'k';
    elseif ~isempty(strfind(fname, 'after'))
        color = 'r';
    else
        color = 'y';
    end
else
    color = C;
end
end

function plotOffset (out, symbols, d, colors, index)
hold on;
j=0;
for i=1:numel(out)
    % set color
    color = getColor(d(i).name, colors{i});
    % set symbol
    j=j+1;
    if j>numel(symbols)
        j=1;
    end
    symbol = symbols{j};
    % set symbol size
    s = getSize(out{i}.final_MI_frac);
    % plot
    plot( 0, out{i}.offsetF(index), 'Marker', symbol, 'Color', color, 'MarkerSize',s);
    plot( 0, out{i}.offsetI(index), 'Marker', '.', 'Color', [1 0 0], 'MarkerSize',8);
end
hold off;
end

function plotRot (out, symbols, d, colors, index)
hold on;
j=0;
for i=1:numel(out)
    % set color
    color = getColor(d(i).name, colors{i});
    % set symbol
    j=j+1;
    if j>numel(symbols)
        j=1;
    end
    symbol = symbols{j};
    % set symbol size
    s = getSize(out{i}.final_MI_frac);
    % plot
    plot( 0, out{i}.rot(index), 'Marker', symbol, 'Color', color, 'MarkerSize',s);
    %plot( 0, out{i}.angle(1), 'Marker', symbols{j}, 'Color', colors{i}, 'MarkerSize',8 );
    plot( 0, out{i}.angle(index), 'Marker', '.', 'Color', [1 0 0], 'MarkerSize',8);
end
hold off;
end




%%
function plot_raw (out, ppath, colors, symbols, d)

%a = numel(out);
f = figure;

subplot(1,3,1);
plotOffset (out, symbols, d, colors, 1);
ylabel('offset in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
plotOffset (out, symbols, d, colors, 2);
ylabel('offset in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
plotOffset (out, symbols, d, colors, 3);
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

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% str = sprintf('Circles are final parameter values. Dots are initial parameter values.');
% text(0.15,0.08,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
% str = sprintf('Black circles are volumes before live sample. Red circles are volumes after live sample.');
% text(0.15,0.04,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_offset.png'];
print(f,str,'-dpng');

% rotation
h = figure;

subplot(1,3,1);
plotRot (out, symbols, d, colors, 1);
ylabel('rotation around dim 1 (radians)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);


subplot(1,3,2);
plotRot (out, symbols, d, colors, 2);
ylabel('rotation around dim 2 (radians)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
plotRot (out, symbols, d, colors, 3);
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

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1);
% str = sprintf('Circles are final parameter values. Dots are initial parameter values.');
% text(0.15,0.08,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
% str = sprintf('Black circles are volumes before live sample. Red circles are volumes after live sample.');
% text(0.15,0.04,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_rotation.png'];
print(h,str,'-dpng');

end

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