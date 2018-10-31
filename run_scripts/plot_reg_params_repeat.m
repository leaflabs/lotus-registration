clear all;
close all;

top_path = '/Users/justin/Desktop/DDLFM/run_1_(distribution_assesment)_on_181010_at_1831';
d = dir([top_path '/**/*.log']);

% Each of 25 data sets was registered 6 times.
% The first 3 times were with different seeds to random number generator.
% The second 3 times were using median registration params from first 3 as coarse registration.
stride = 6;

%% Gather data from files
tmp_f = 'tmp.mat';
fname = [top_path '/' tmp_f];
if ~exist(fname,'file')
    offset = {};
    rot = {};
    offsetmc = {};
    rotmc = {};
    % for each RNG seed, i
    for i=1:3
        % capture registration parameters using default coarse registration
        [off, r, centroid] = read_set( d, i:stride:numel(d) );
        offset{i} = off;
        rot{i} = r;
        % capture registration parameters using median params for coarse registration
        [off , r] = read_set( d, i+3:stride:numel(d) );
        offsetmc{i} = off;
        rotmc{i} = r;
    end
    save(fname,'offset','rot','offsetmc','rotmc','centroid');
else
    load(fname);
end

%ground truth
offset_gt = [-5 -10 15];
rot_gt = [-87 2 -4]*pi/180;

%% main
fname = [top_path '/reg_params_across_seeds_default_coarse.png'];
plot_seeds (offset, rot, 'default coarse registration', fname);
fname = [top_path '/reg_params_across_seeds_median_coarse.png'];
plot_seeds (offsetmc, rotmc, 'coarse registration using median params', fname);

% calculate registration parameter error
% for each RNG seed, i
offset_err = {};
offsetmc_err = {};
rot_err = {};
rot_mc_err = {};
for i=1:3
    offset_err{i} = offset{i}-offset_gt;
    offsetmc_err{i} = offsetmc{i}-offset_gt;
    rot_err{i} = rot{i}-rot_gt;
    rotmc_err{i} = rotmc{i}-rot_gt;
end

print_tables (offset_err, offsetmc_err, rot_err, rotmc_err, top_path, centroid);

% for each RNG seed, i
for i=1:3
    fname = [top_path sprintf('/offset_rot_error_%d.png',i)];
    plot_error( offset{i}, offsetmc{i}, rot{i}, rotmc{i}, ...
    offset_gt, rot_gt, fname, i );
end

%%
function print_tables (offset, offsetmc, rot, rotmc, top_path, centroid)
fprintf('\n%6s%25s%25s%25s%25s%25s%25s\n','seed','off_err_dim1[um]','off_err_dim2[um]','off_err_dim3[um]','rot_err_dim1[deg]','rot_err_dim2[deg]','rot_err_dim3[deg]');
tabulate (offset, offsetmc, rot, rotmc);
r = sqrt(sum(centroid.*centroid));
fprintf('\nWorst-case moment arm = %3.1f um\n',r);
rot_um = {};
for i=1:3
    rot_um{i} = rot{i}*pi/180*r;
    rotmc_um{i} = rotmc{i}*pi/180*r;
end
fprintf('\n%6s%25s%25s%25s%25s%25s%25s\n','seed','off_err_dim1[um]','off_err_dim2[um]','off_err_dim3[um]','rot_err_dim1[um]','rot_err_dim2[um]','rot_err_dim3[um]');
tabulate (offset, offsetmc, rot_um, rotmc_um);
end


function tabulate (offset, offsetmc, rot, rotmc)
corr0 = '';
corr1 = '';
for seed=1:3
    msg0 = sprintf('%6s',num2str(seed));
    msg1 = msg0;
    for dim=1:3
        [str0, str1] = get_str( offset{seed}(:,dim), offsetmc{seed}(:,dim) );
        msg0 = [msg0 sprintf('%25s',str0)];
        msg1 = [msg1 sprintf('%25s',str1)];
    end
    for dim=1:3
        [str0, str1] = get_str( rot{seed}(:,dim)*180/pi, rotmc{seed}(:,dim)*180/pi );
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


function plot_error ( offset, offset_mc, rot, rot_mc, offset_gt, rot_gt, fname, seed )

offset_err = offset-offset_gt;
offset_mc_err = offset_mc-offset_gt;
rot_err = rot-rot_gt;
rot_mc_err = rot_mc-rot_gt;

a = size(offset_err);
f = figure;

subplot(1,2,1);
hold on;
plot(1*ones(1,a(1)),offset_err(:,1),'ko','LineWidth',2);
plot(2*ones(1,a(1)),offset_err(:,2),'ko','LineWidth',2);
plot(3*ones(1,a(1)),offset_err(:,3),'ko','LineWidth',2);
plot(1*ones(1,a(1)),offset_mc_err(:,1),'ro','LineWidth',2);
plot(2*ones(1,a(1)),offset_mc_err(:,2),'ro','LineWidth',2);
plot(3*ones(1,a(1)),offset_mc_err(:,3),'ro','LineWidth',2);
hold off;
xlabel('translation dims');
ylabel('error [um]');
ylim([-20 20]);

subplot(1,2,2);
hold on;
plot(4*ones(1,a(1)),rot_err(:,1)*180/pi,'ko','LineWidth',2);
plot(5*ones(1,a(1)),rot_err(:,2)*180/pi,'ko','LineWidth',2);
plot(6*ones(1,a(1)),rot_err(:,3)*180/pi,'ko','LineWidth',2);
plot(4*ones(1,a(1)),rot_mc_err(:,1)*180/pi,'ro','LineWidth',2);
plot(5*ones(1,a(1)),rot_mc_err(:,2)*180/pi,'ro','LineWidth',2);
plot(6*ones(1,a(1)),rot_mc_err(:,3)*180/pi,'ro','LineWidth',2);
hold off;
xlabel('rotation dims');
ylabel('error [degrees]');
title(sprintf('random number generator seed = %d',seed));

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
