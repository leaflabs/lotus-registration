function param = read_logs (param)
%% read log files
pathNpattern = sprintf('%s*FrameNumber_%04d_*.log',param.inter,param.m);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end
content = fileread([logs(1).folder '/' logs(1).name]);

eval(['param.transM = ' get_match('trans',content) ';']);
eval(['param.rotM = ' get_match('rot',content) ';']);
eval(['param.centroidM = ' get_match('centroid',content) ';']);
eval(['param.clipM = ' get_match('clip',content) ';']);

pathNpattern = sprintf('%s*FrameNumber_%04d_*.log',param.inter,param.n);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end
content = fileread([logs(1).folder '/' logs(1).name]);

eval(['param.transN = ' get_match('trans',content) ';']);
eval(['param.rotN = ' get_match('rot',content) ';']);
eval(['param.centroidN = ' get_match('centroid',content) ';']);
eval(['param.clipN = ' get_match('clip',content) ';']);

if ~isempty(find( (param.clipN-param.clipM) ~= 0))
    disp('WTF?');
    keyboard
end
param.clip = param.clipN; % equal to clipM
if ~isempty(find( (param.centroidN-param.centroidM) ~= 0))
    disp('WTF?');
    keyboard
end
param.centroid = param.centroidN ...
        + [param.clip(1)*param.voxel_y param.clip(3)*param.voxel_x 0];
end

function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{1}{1};
end


