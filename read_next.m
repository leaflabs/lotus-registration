function next = read_next (ddir, prefix, suffix)
%% read log files
pathNpattern = sprintf('%s%s*%s',ddir,prefix,suffix);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end

v = 0;
for i=1:numel(logs)
    t = logs(i).name(1:end-4);
    j = str2num(t(end-3:end));
    if j>v
        v = j;
    end
end
next = v+1;
