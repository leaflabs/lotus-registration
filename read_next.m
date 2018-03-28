function next = read_next (ddir, prefix, suffix)
%% read log files
pathNpattern = sprintf('%s%s*%s',ddir,prefix,suffix);
logs = dir(pathNpattern);
if isempty(logs)
    fprintf('Warning. No matches found for %s\n',pathNpattern);
    if ~exist(ddir,'dir')
        fprintf('Folder does not exist:%s\n',ddir);
        exit;
    else
        % folder found so assume fresh start
        next = 1;
    end
else
    % process matches
    v = 0;
    for i=1:numel(logs)
        t = logs(i).name(1:end-4);
        j = str2num(t(end-3:end));
        if j>v
            v = j;
        end
    end
    next = v+1;
end
