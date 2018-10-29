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
    % RNG seed 1
    [offset1 , rot1] = read_set( d, 1:stride:numel(d) );
    [offset1_mc , rot1_mc] = read_set( d, 4:stride:numel(d) );
    % RNG seed 2
    [offset2 , rot2] = read_set( d, 2:stride:numel(d) );
    [offset2_mc , rot2_mc] = read_set( d, 5:stride:numel(d) );
    % RNG seed 3
    [offset3 , rot3] = read_set( d, 3:stride:numel(d) );
    [offset3_mc , rot3_mc] = read_set( d, 6:stride:numel(d) );
    save(fname);
else
    load(fname);
end

offset_ground_truth = [-5 -10 15];
rot_ground_truth = [-87 2 -4]*pi/180;

plot_error( offset1, offset1_mc, rot1, rot1_mc, ...
    offset_ground_truth, rot_ground_truth, ...
    [top_path '/offset_rot_error_1.png'] );

plot_error( offset2, offset2_mc, rot2, rot2_mc, ...
    offset_ground_truth, rot_ground_truth, ...
    [top_path '/offset_rot_error_2.png'] );

plot_error( offset3, offset3_mc, rot3, rot3_mc, ...
    offset_ground_truth, rot_ground_truth, ...
    [top_path '/offset_rot_error_3.png'] );

%%

function plot_error ( offset, offset_mc, rot, rot_mc, offset_gt, rot_gt, fname )

offset_err = offset-offset_gt;
offset_mc_err = offset_mc-offset_gt;
rot_err = rot-rot_gt;
rot_mc_err = rot_mc-rot_gt;

a = size(offset_err);
f = figure;

subplot(1,2,1);
hold on;
plot(1*ones(1,a(1)),offset_err(:,1),'ko');
plot(2*ones(1,a(1)),offset_err(:,2),'ko');
plot(3*ones(1,a(1)),offset_err(:,3),'ko');
plot(1*ones(1,a(1)),offset_mc_err(:,1),'ro');
plot(2*ones(1,a(1)),offset_mc_err(:,2),'ro');
plot(3*ones(1,a(1)),offset_mc_err(:,3),'ro');
hold off;
xlabel('translation dims');
ylabel('error [um]');
ylim([-20 20]);

subplot(1,2,2);
hold on;
plot(4*ones(1,a(1)),rot_err(:,1)*180/pi,'ko');
plot(5*ones(1,a(1)),rot_err(:,2)*180/pi,'ko');
plot(6*ones(1,a(1)),rot_err(:,3)*180/pi,'ko');
plot(4*ones(1,a(1)),rot_mc_err(:,1)*180/pi,'ro');
plot(5*ones(1,a(1)),rot_mc_err(:,2)*180/pi,'ro');
plot(6*ones(1,a(1)),rot_mc_err(:,3)*180/pi,'ro');
hold off;
xlabel('rotation dims');
ylabel('error [degrees]');

print(f,fname,'-dpng');
end



function [offset , rot] = read_set ( d, vec )
offset = [];
rot = [];
for i=vec
    fpath = [d(i).folder '/' d(i).name];
    fprintf('Reading log file: %s\n',fpath);
    out = read_log (fpath);
    offset = [offset ; out.offsetF];
    rot = [rot ; out.rot];
end
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
