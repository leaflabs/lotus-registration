clear all;
close all;

top_path = '/Users/justin/Desktop/DDLFM/run_1_(distribution_assesment)_on_181010_at_1831';
d = dir([top_path '/**/*.log']);

% Each of 25 data sets was registered 6 times.
% The first 3 times were with different seeds to random number generator.
% The second 3 times were using median registration params from first 3 as coarse registration.
stride = 2;
seeds = 1;

%% Gather data from files
tmp_f = 'tmp.mat';
fname = [top_path '/' tmp_f];
if ~exist(fname,'file')
    offset = {};
    rot = {};
    offsetmc = {};
    rotmc = {};
    % for each RNG seed, i
    for i=1:seeds
        % capture registration parameters using default coarse registration
        [off, r, centroid] = read_set( d, i:stride:numel(d) );
        offset{i} = off;
        rot{i} = r;
        % capture registration parameters using median params for coarse registration
        [off , r] = read_set( d, i+seeds:stride:numel(d) );
        offsetmc{i} = off;
        rotmc{i} = r;
    end
    [offset_median , rot_median] = get_median( d, 1:2:numel(d) );
    [offsetmc_median , rotmc_median] = get_median( d, 2:2:numel(d) );
    save(fname,'offset','rot','offsetmc','rotmc','centroid','offset_median','rot_median','offsetmc_median','rotmc_median');
else
    load(fname);
end

%% main

%ground truth
%offset_gt = [-5 -10 15];
offset_gt = [14 16 5]; % um
g = sprintf('%d ', offset_gt);
offset_gt_str = sprintf('Assume gt = [%s] um', g); 
%rot_gt = [87 2 -4]*pi/180;
rot_gt_deg = [87 -4 -2]; % radians
rot_gt_rad = rot_gt_deg * pi/180; % radians
g = sprintf('%d ', rot_gt_deg);
rot_gt_str = sprintf('Assume gt = [%s] deg', g);



if seeds>1
    fname = [top_path '/reg_params_across_seeds_default_coarse.png'];
    plot_seeds (offset, rot, 'default coarse registration', fname);
    fname = [top_path '/reg_params_across_seeds_median_coarse.png'];
    plot_seeds (offsetmc, rotmc, 'use median param values for coarse registration', fname);
else
    fname = [top_path '/reg_params_across_data_sets_default_coarse.png'];
    plot_raw (offset{1}, rot{1}, 'default coarse registration', fname);
    fname = [top_path '/reg_params_across_data_sets_median_coarse.png'];
    plot_raw (offsetmc{1}, rotmc{1}, 'use median param values for coarse registration', fname);
end

 
% calculate registration parameter error
% for each RNG seed, i
offset_err = {};
offsetmc_err = {};
rot_err_rad = {};
rotmc_err_rad = {};
rot_err_um = {};
rotmc_err_um = {};
r = sqrt(sum(centroid.*centroid));
fprintf('\nWorst-case moment arm = %3.1f um\n',r);
for i=1:seeds
    offset_err{i}    = offset{i}   - offset_gt;
    offsetmc_err{i}  = offsetmc{i} - offset_gt;
    rot_err_rad{i}   = rot{i}      - rot_gt_rad;
    rotmc_err_rad{i} = rotmc{i}    - rot_gt_rad;
    rot_err_um{i}    = rot_err_rad{i}   * r;
    rotmc_err_um{i}  = rotmc_err_rad{i} * r;
end

print_tables (offset_err, offsetmc_err, rot_err_rad, rotmc_err_rad, rot_err_um, rotmc_err_um,seeds);

% for each RNG seed, i
for i=1:seeds
    fname = [top_path sprintf('/offset_rot_error_%d.png',i)];
    plot_error( offset_err{i}, offsetmc_err{i}, rot_err_rad{i}, rotmc_err_rad{i}, ...
    rot_err_um{i}, rotmc_err_um{i}, r, fname, i, offset_gt_str, rot_gt_str );
end


% centroids
rnd1_offset_err = offset_median - offset_gt;
rnd2_offset_err = offsetmc_median - offset_gt;
rnd1_rot_err_rad = rot_median - rot_gt_rad;
rnd2_rot_err_rad = rotmc_median - rot_gt_rad;
rnd1_rot_err_um = rnd1_rot_err_rad * r;
rnd2_rot_err_um = rnd2_rot_err_rad * r;

fprintf('\n\n');
g = sprintf('%3.1f ', offset_gt);
h = sprintf('%3.1f ', rot_gt_rad*180/pi);
fprintf('                         Ground truth: offset = [%s] um,        rotation = [%s] deg\n', g, h);
g = sprintf('%3.1f ', offset_median);
h = sprintf('%3.1f ', rot_median*180/pi);
fprintf('   Default coarse registration: median offset = [%s] um, median rotation = [%s] deg\n', g, h);
g = sprintf('%3.1f ', rnd1_offset_err);
h = sprintf('%3.1f ', rnd1_rot_err_rad*180/pi);
i = sprintf('%3.1f ', rnd1_rot_err_um);
fprintf('                                        error = [%s] um,           error = [%s] deg ~= [%s] um\n', g, h, i);
g = sprintf('%3.1f ', offsetmc_median);
h = sprintf('%3.1f ', rotmc_median*180/pi);
fprintf('Median for coarse registration: median offset = [%s] um, median rotation = [%s] deg\n', g, h);
g = sprintf('%3.1f ', rnd2_offset_err);
h = sprintf('%3.1f ', rnd2_rot_err_rad*180/pi);
i = sprintf('%3.1f ', rnd2_rot_err_um);
fprintf('                                        error = [%s] um,           error = [%s] deg ~= [%s] um\n', g, h, i);
fprintf('\n\n');





%%
function [median_offset , median_rot] = get_median ( d, vec )
offset = [];
rot = [];
for i=vec
    fpath = [d(i).folder '/' d(i).name];
    fprintf('Reading log file: %s\n',fpath);
    out = read_log (fpath);
    offset = [offset ; out.offsetF];
    rot = [rot ; out.rot];
end
median_offset = median(offset);
median_rot = median(rot);
end


function print_tables (offset_err, offsetmc_err, rot_err_rad, rotmc_err_rad, rot_err_um, rotmc_err_um, seeds)
fprintf('\n%6s%25s%25s%25s%25s%25s%25s\n','seed',...
    'off_err_dim1[um]','off_err_dim2[um]','off_err_dim3[um]',...
    'rot_err_dim1[deg]','rot_err_dim2[deg]','rot_err_dim3[deg]');
tabulate (offset_err, offsetmc_err, rot_err_rad, rotmc_err_rad, seeds);
fprintf('\n%6s%25s%25s%25s%25s%25s%25s\n','seed',...
    'off_err_dim1[um]','off_err_dim2[um]','off_err_dim3[um]',...
    'rot_err_dim1[um]','rot_err_dim2[um]','rot_err_dim3[um]');
tabulate (offset_err, offsetmc_err, rot_err_um, rotmc_err_um, seeds);
end


function tabulate (offset_err, offsetmc_err, rot_err_rad, rotmc_err_rad, seeds)
corr0 = '';
corr1 = '';
for seed=1:seeds
    msg0 = sprintf('%6s',num2str(seed));
    msg1 = msg0;
    for dim=1:3
        [str0, str1] = get_str( offset_err{seed}(:,dim), offsetmc_err{seed}(:,dim) );
        msg0 = [msg0 sprintf('%25s',str0)];
        msg1 = [msg1 sprintf('%25s',str1)];
    end
    for dim=1:3
        [str0, str1] = get_str( rot_err_rad{seed}(:,dim)*180/pi, rotmc_err_rad{seed}(:,dim)*180/pi );
        msg0 = [msg0 sprintf('%25s',str0)];
        msg1 = [msg1 sprintf('%25s',str1)];
    end
    corr0 = [corr0 sprintf('%s\n',msg0)];
    corr1 = [corr1 sprintf('%s\n',msg1)];
end
fprintf('%s',corr0);
fprintf('\n');
fprintf('%s',corr1);
end

function [str0,str1] = get_str ( val0, val1 )
pd = fitdist(val0,'Normal');
pdmc = fitdist(val1,'Normal');
str0 = sprintf('mu: %3.1f -> %3.1f',pd.mu, pdmc.mu);
str1 = sprintf('sigma: %3.1f -> %3.1f',pd.sigma, pdmc.sigma);
end

function plot_seeds (offset, rot, mytitle, fname)
%% plot variation across seeds
n_offset = [offset{1}; offset{2}; offset{3}];
n_rot = [rot{1}; rot{2}; rot{3}];
% derive normalization value for each of 6 registration parameters
offset_m = median(n_offset);
rot_m = median([rot{1}; rot{2}; rot{3}]);
% normalize all data
n_offset = n_offset./offset_m;
n_rot = n_rot./rot_m;

f = figure;

a = size(offset{1});
valmat = [];
index = [];

b = 1;
c = a(1);
vals = [ n_offset(b:c,1); n_offset(b:c,2); n_offset(b:c,3); ...
    n_rot(b:c,1); n_rot(b:c,2); n_rot(b:c,3)];
valmat = [valmat vals];
index = [index 1*ones(numel(vals),1)];

b = 1+a(1);
c = 2*a(1);
vals = [ n_offset(b:c,1); n_offset(b:c,2); n_offset(b:c,3); ...
    n_rot(b:c,1); n_rot(b:c,2); n_rot(b:c,3)];
valmat = [valmat vals];
index = [index 2*ones(numel(vals),1)];

b = 1+2*a(1);
c = 3*a(1);
vals = [ n_offset(b:c,1); n_offset(b:c,2); n_offset(b:c,3); ...
    n_rot(b:c,1); n_rot(b:c,2); n_rot(b:c,3)];
valmat = [valmat vals];
index = [index 3*ones(numel(vals),1)];

plot(index',valmat','ko-','LineWidth',2);
xlabel('RNG seeds');
ylabel('registration parameters normalized by median');
xlim([0 4]);
title(mytitle);
print(f,fname,'-dpng');
end

function plot_raw (offset, rot, mytitle, fname)

a = size(offset);

f = figure;
subplot(1,2,1);
hold on;
plot(1*ones(1,a(1)),offset(:,1),'ko','LineWidth',2);
plot(2*ones(1,a(1)),offset(:,2),'ko','LineWidth',2);
plot(3*ones(1,a(1)),offset(:,3),'ko','LineWidth',2);
hold off;
xlabel('translation dims');
ylabel('offset [um]');
xlim([0 4]);

subplot(1,2,2);
hold on;
plot(1*ones(1,a(1)),rot(:,1)*180/pi,'ko','LineWidth',2);
plot(2*ones(1,a(1)),rot(:,2)*180/pi,'ko','LineWidth',2);
plot(3*ones(1,a(1)),rot(:,3)*180/pi,'ko','LineWidth',2);
hold off;
xlabel('rotation dims');
ylabel('angle [degrees]');
xlim([0 4]);

title(mytitle);

print(f,fname,'-dpng');
end


function plot_error ( offset_err, offsetmc_err, ...
    rot_err_rad, rotmc_err_rad, ...
    rot_err_um, rotmc_err_um, ...
    r, fname, seed, offset_gt_str, rot_gt_str )

a = size(offset_err);
f = figure;
x = [0 4];

subplot(1,3,1);
hold on; 
plot(1*ones(1,a(1)),offsetmc_err(:,1),'ko','LineWidth',2);
plot(2*ones(1,a(1)),offsetmc_err(:,2),'ko','LineWidth',2);
plot(3*ones(1,a(1)),offsetmc_err(:,3),'ko','LineWidth',2);
% y = ylim;
% plot(1*ones(1,a(1)),offset_err(:,1),'ko','LineWidth',2);
% plot(2*ones(1,a(1)),offset_err(:,2),'ko','LineWidth',2);
% plot(3*ones(1,a(1)),offset_err(:,3),'ko','LineWidth',2);
hold off;
xlabel('translation dims');
ylabel('error [um]');
%ylim(y);
xlim(x);
title(offset_gt_str);

subplot(1,3,2);
hold on;
plot(1*ones(1,a(1)),rotmc_err_rad(:,1)*180/pi,'ko','LineWidth',2);
plot(2*ones(1,a(1)),rotmc_err_rad(:,2)*180/pi,'ko','LineWidth',2);
plot(3*ones(1,a(1)),rotmc_err_rad(:,3)*180/pi,'ko','LineWidth',2);
% y = ylim;
% plot(1*ones(1,a(1)),rot_err_rad(:,1)*180/pi,'ko','LineWidth',2);
% plot(2*ones(1,a(1)),rot_err_rad(:,2)*180/pi,'ko','LineWidth',2);
% plot(3*ones(1,a(1)),rot_err_rad(:,3)*180/pi,'ko','LineWidth',2);
hold off;
xlabel('rotation dims');
ylabel('error [degrees]');
%title(sprintf('random number generator seed = %d',seed));
%ylim(y);
xlim(x);
title(rot_gt_str);

subplot(1,3,3);
hold on;
plot(1*ones(1,a(1)),rotmc_err_um(:,1),'ko','LineWidth',2);
plot(2*ones(1,a(1)),rotmc_err_um(:,2),'ko','LineWidth',2);
plot(3*ones(1,a(1)),rotmc_err_um(:,3),'ko','LineWidth',2);
% y = ylim;
% plot(1*ones(1,a(1)),rot_err_um(:,1),'ko','LineWidth',2);
% plot(2*ones(1,a(1)),rot_err_um(:,2),'ko','LineWidth',2);
% plot(3*ones(1,a(1)),rot_err_um(:,3),'ko','LineWidth',2);
hold off;
xlabel('rotation dims');
ylabel('error [um]');
%ylim(y);
xlim(x);
title(sprintf('Worst-case moment arm = %4.0f um',r));


print(f,fname,'-dpng');
end



function [offset , rot, centroid] = read_set ( d, vec )
offset = [];
rot = [];
for i=vec
    fpath = [d(i).folder '/' d(i).name];
    fprintf('Reading log file: %s\n',fpath);
    out = read_log (fpath);
    offset = [offset ; out.offsetF];
    rot = [rot ; out.rot];
end
centroid = out.centroid;
end


function out = read_log (f)
out = [];
content = fileread(f);
eval(['out.trans = ' get_match('trans',content) ';']);
%eval(['out.angle = ' get_match('angle',content) ';']);
eval(['out.centroid = ' get_match('centroid',content) ';']);
%eval(['out.offsetI = ' get_match('offset',content) ';']);
out.offsetF = out.trans - out.centroid;
eval(['out.rot = ' get_match('rot',content) ';']);
end


function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{end}{1};
end
