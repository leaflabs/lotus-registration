clear all;
close all;

ppath = '/home/jkinney/Desktop/DDLFM/fish/20180123/interpolate';
pathNpattern = sprintf('%s/*.log',ppath);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end

r1 = [];
r2 = [];
r3 = [];
n = numel(logs);
for i=1:n
    fprintf('%d of %d\n',i,n);
    f = [logs(i).folder '/' logs(i).name]
    rot = read_log (f);
    r1 = [r1;rot(:,1)'];
    r2 = [r2;rot(:,2)'];
    r3 = [r3;rot(:,3)'];
end

h = figure('units','normalized','outerposition',[0 0 1 1]);

hold on;
for i=1:n
    fprintf('%d of %d\n',i,n);
    plot(r1(i,:));
end
hold off;
keyboard

% subplot(1,3,1);
% plot(arc(:,1));
% xlabel('iteration');
% ylabel('arc about dim 1 (um)');
% 
% subplot(1,3,2);
% plot(arc(:,1));
% ylabel('arc about dim 2 (um)');
% xlabel('iteration');
% 
% subplot(1,3,3);
% plot(arc(:,1));
% ylabel('arc about dim 3 (um)');
% xlabel('iteration');


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'Arc amplitude due to rotation during simulated annealing';
text(0.25,0.97,str,'FontSize',12,'Color',[0 0 0],'Interpreter','none');
text(0.25,0.94,ppath,'FontSize',12,'Color',[0 0 0],'Interpreter','none');


str = [ppath '/arc_amplitude.png'];
print(h,str,'-dpng');

function rot = read_log (f)
%% read log file
out = [];

content = fileread(f);
C = strsplit(content,'\n');
% fid = fopen(f);
% content = textscan(fid,'%s\n');
% fclose(fid);

expr = ['.* r = [^\n]*trans[^\n]*'];
rot = get_rot_match(expr,C);
%eval(['xlim_1d = ' get_match3('xlim_1d',content) ';']);
%eval(['voxel_x = ' get_match1('voxel_x',content) ';']);
%eval(['voxel_y = ' get_match1('voxel_y',content) ';']);

%half_span = 0.5 * [ xlim_1d(1)*voxel_y   xlim_1d(2)*voxel_x ];
%radius = sqrt( sum( half_span .* half_span ) );

%arc = radius * rot;

end



function out = get_rot_match (expr, C)
out = [];
for i=1:numel(C)
    line = C{i};
    if ~isempty(line)
        match = regexp(line,expr,'match');
        if ~isempty(match)
            A = strsplit(match{1},'test MI');
            B = strsplit(A{1},',');
            eval([B{end} ';']);
            out = [out;rot];
        end
    end
end
end

function out = get_match3 (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{1}{1};
end


function out = get_match1 (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = strsplit(match{1},':');
out = s{2};
end
