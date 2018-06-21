clear all;
close all;

voxel = 1; % um
volume_size = 150; % voxels

[~,hostname] = system('hostname')

if ~isempty(strfind(hostname, 'Justins-Mac'))
	ppath = '/Users/justin/Desktop/DDLFM';
elseif ~isempty(strfind(hostname, 'willis'))
	ppath = '/home/jkinney/Desktop/DDLFM';
else
	disp('WTF?');
    keyboard
end

%ddlfm_path = 'bead/20180328/registration/';
ddlfm_path = 'Transformation-Phantoms-master/analysis/registration/repeat_90_90_0_D/';
opath = [ppath '/' ddlfm_path];

savePath = [opath 'scan_results/'];
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
fname = sprintf('%sscan_corbans_phantom_%s.log',savePath, timestamp);
diary(fname)
tic

tmp_f = 'tmp.mat';

%% Gather data from files
if ~exist([savePath tmp_f],'file')
    % % row col z
    % % y x z
    % % dim 1 2 3
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

colcell = {'target_offset'
    'initial_offset'
    'estimated_offset'
    'final_offset'
    'offset_error_um'
    'offset_error_voxel'
    'target_rot'
    'initial_rot'
    'final_rot'
    'rot_error_deg'
    'mean_rot_error_voxel'};

padsize = 25;
header = pad('index',7);
for i=1:numel(colcell)
    header = [header pad(colcell{i},padsize)];
end
fprintf('offset in units of micrometers. Angle in units of degrees.\n');
fprintf('%s\n',header);

mr = estimate_mean_radius(volume_size); % voxels

for i=1:numel(out)
    out{i}.offset_error_um = out{i}.final_offset-out{i}.target_offset;
    out{i}.final_rot = out{i}.final_rot * 180/pi;
    out{i}.initial_rot = out{i}.initial_rot * 180/pi;
    out{i}.rot_error_deg = out{i}.final_rot-out{i}.target_rot;
    out{i}.offset_error_voxel = out{i}.offset_error_um * voxel;
    out{i}.mean_rot_error_voxel = mr * out{i}.rot_error_deg * pi/180;
    row = pad(num2str(i),7);
    for j=1:numel(colcell)
        eval(['tmp = out{i}.' colcell{j} ';'])
        mystr = sprintf('%2.1f  %2.1f  %2.1f',tmp(1),tmp(2),tmp(3));
        row = [row pad(mystr,padsize)];
    end
    fprintf('%s\n',row);
end

diary off;

%%
function out = estimate_mean_radius (volume_size)
N = round(volume_size/2);
sum = 0;
for r = 1:N
    pop = (2*r+1)^3 - (2*r-1)^3;
    sum = sum + r*pop;
end
out = sum / volume_size^3;
end


%%
function out = read_log (f)
out = [];
content = fileread(f);

% offset
target_offset = getOffsetFromName(f);
eval(['out.target_offset = [' target_offset '];']);
trans_matches =  get_match(' trans:',content);
eval(['out.initial_offset = [' trans_matches{1}{1} '];']);
estimated_matches =  get_match('estimated offsets = ',content);
eval(['out.estimated_offset = ' estimated_matches{1}{1} ';']);
centroid_matches =  get_match(' centroid:',content);
eval(['centroid = [' centroid_matches{2}{1} '];']);
eval(['final_trans = [' trans_matches{2}{1} '];']);
out.final_offset = final_trans - centroid;

% rot
target_rot = getRotFromName(f);
eval(['out.target_rot = [' target_rot '];']);
angle_matches =  get_match(' angle:',content);
eval(['out.initial_rot = [' angle_matches{1}{1} '];']);
rot_matches =  get_match(' rot:',content);
eval(['out.final_rot = [' rot_matches{2}{1} '];']);

end

%%
function out = get_match (str, content)
expr = ['[^\n]*' str '[^\n]*'];
match = regexp(content,expr,'match');
out = regexp(match,'\[.*\]','match');
end

%%
function out = getRotFromName (log_filename)
content = log_filename;
expr = 'rot_[^\n]*_trans';
match = regexp(content,expr,'match');
% Assume example: Recon3D_after_10_Mono
vector = match{1}(5:end-6);
C = strsplit(vector,'_');
out = [C{1} ' ' C{2} ' ' C{3}];
end

function out = getOffsetFromName (log_filename)
content = log_filename;
expr = 'trans_[^\n]*__';
match = regexp(content,expr,'match');
% Assume example: Recon3D_after_10_Mono
vector = match{1}(7:end-2);
C = strsplit(vector,'_');
out = [C{1} ' ' C{2} ' ' C{3}];
end




