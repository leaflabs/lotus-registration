clear all;
close all;

%% load parameters
param = init_param();

%%
[~,hostname] = system('hostname')
param.hostname = strtrim(hostname);

setup_par

if ~isempty(strfind(param.hostname, 'Justins-Mac'))
	param.ppath = '/Users/justin/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180328'];
	param.opath = param.ipath;
    param.parallel = false;
elseif ~isempty(strfind(param.hostname, 'willis'))
	param.ppath = '/home/jkinney/Desktop/DDLFM';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = [param.ppath '/bead/20180328'];
	param.opath = param.ipath;
    param.parallel = false;
else
	param.ppath = '/home/jkinney';
	addpath([param.ppath '/lotus-registration']);
	param.ipath = '/om/user/jkinney/DLFM/bead/20180328';
	param.opath = param.ipath;
    param.parallel = false;
end
param.savePath = [param.opath '/registration_20180520_170824/analysis/'];

if ~exist(param.savePath,'dir')
    status = mkdir(param.savePath);
    if status == 1
        disp(['Created folder: ' param.savePath]);
    else
        disp(['Error attempting to create folder:' param.savePath]);
        status
        exit;
    end
end

%%
timestamp = [datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%sshuffle_bead_before_after_%s.log',param.savePath, timestamp);
diary(fname)
tic

%%
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0;     % um
param.justCalcMI = true;

%% Gather registration parameters from files
if 0>1
    % % row col z
    % % y x z
    % % dim 1 2 3
    beads = {};
    out = {};
    p = [param.ipath '/registration_20180520_170824/*.log'];
    d = dir(p);
    fprintf('%d log files found for\n%s\n\n',numel(d),p);
    for i=1:length(d)
        fname = [d(i).folder '/' d(i).name]
        % read log to get registration parameters for volume
        %lf = [fname(1:end-18) '.log'];
        out{i} = read_log (fname);
        out{i}
    end
    save([param.savePath 'tmp_shuffle.mat']);
else
    load([param.savePath 'tmp_shuffle.mat'],'out','d');
end

param.inputFilePath1 = [param.ipath '/horizontal/Reconstructed/'];
param.inputFilePath2 = [param.ipath '/vertical/Reconstructed/'];
param.rapid = 0;
param.justCalcMI = true;

%% apply registration parameters
N = numel(out);

if 0>1
    MIvec = zeros(N,N);
    for i=1:N
        % apply registration parameter set i
        param.rot = out{i}.rot;
        %     param.trans = out{i}.trans + [ out{i}.clip(1)*param.voxel_y  ...
        %         out{i}.clip(3)*param.voxel_x   out{i}.clip(5)*param.voxel_z ];
        param.trans = out{i}.trans;
        % to all N volume pairs
        for j=1:N
            fprintf('\n####\n\nRegistration parameter set %d of %d\n',i,N);
            fprintf('Volume pair %d of %d\n',j,N);
            fprintf('Derived from data set %s\n', getInputFileName( d(i).name ));
            fprintf('Apply to data set %s\n',getInputFileName( d(j).name ));
            fprintf('rot = [%f %f %f]\n',param.rot(1),param.rot(2),param.rot(3));
            fprintf('trans = [%f %f %f]\n',param.trans(1),param.trans(2),param.trans(3));
            param.inputFileName = {getInputFileName( d(j).name )};
            MI = register(param);
            fprintf('MI = %f\n',MI);
            close all;
            MIvec(i,j) = MI;
        end
    end
    save([param.savePath 'tmp_shuffle_MI.mat']);
else
    load([param.savePath 'tmp_shuffle_MI.mat'],'MIvec');
end

%% Plot

f = figure;
hold on;
for j=1:N
    for i=1:N
        if i==j
            plot(j,MIvec(i,j),'ro');
        else
            plot(j,MIvec(i,j),'bo');    
        end
    end
end
hold off;
ylabel('bead position in dim 1 (um)');
set(gca,'xtick',[]);
ylabel('MI');
str = 'shuffle_bead_MI.png';
title({param.savePath;str;'red = match, blue = shuffle'},'Interpreter','None');
%subtitle('red = match, black = shuffle');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
for i=1:N
    name = getDataSetName(d(i).name);
    text(0.165+(i-1)*0.039,0.01,name,'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
    %text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
end

str = [param.savePath 'shuffle_bead_MI.png'];
print(f,str,'-dpng');


f = figure;
hold on;
% for each volume pair
for j=1:N
    % for each set of registration parameters
    for i=1:N
        if i==j
            plot(j,MIvec(i,j)/MIvec(j,j),'ro');
        else
            plot(j,MIvec(i,j)/MIvec(j,j),'bo');    
        end
    end
end
hold off;
ylabel('bead position in dim 1 (um)');
set(gca,'xtick',[]);
ylabel('Normalized MI');
str = 'shuffle_bead_MI_normalized.png';
title({param.savePath;str;'red = match, blue = shuffle'},'Interpreter','None');
%subtitle('red = match, black = shuffle');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
for i=1:N
    name = getDataSetName(d(i).name);
    text(0.165+(i-1)*0.039,0.01,name,'FontSize',9,'Color',[0 0 0],'Interpreter','none','Rotation',90);
    %text(0.3,0.9-(i-1)*0.04,sets{i},'FontSize',9,'Color',colors{i},'Interpreter','none');
end

str = [param.savePath 'shuffle_bead_MI_normalized.png'];
print(f,str,'-dpng');


%% Functions

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

function title = plotOffset (out, symbols, d, colors, index)
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
title = 'before = BLACK, after = RED';
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
mytitle = plotOffset (out, symbols, d, colors, 2);
ylabel('offset in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end
title(mytitle);

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
title(mytitle);

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
eval(['out.clip = ' get_match('clip',content) ';']);
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