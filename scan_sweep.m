clear all;
close all;

%ppath = '/Users/justin/Desktop/DLFM/fish/20171221/repeat';
ppath = '/home/jkinney/Desktop/DLFM/fish/20171221/sweep_rot';
pathNpattern = sprintf('%s/*.log',ppath);
logs = dir(pathNpattern);
if isempty(logs)
    disp('WTF?');
    fprintf('%s\n',pathNpattern);
    keyboard
end

oI = [];
oF = [];
rI = [];
rF = [];
n = numel(logs);
for i=1:n
    fprintf('%d of %d\n',i,n);
    f = [logs(i).folder '/' logs(i).name];
    out = read_log (f);
    oI = [oI;out.offsetI];
    oF = [oF;out.offsetF];
    rF = [rF;out.rot];
    rI = [rI;out.angle];
end

plot_raw(oI, oF, rI, rF, out, ppath);
%plot_histo(o, r, ppath);




%%
function plot_histo (o,r, ppath)

a = size(o);
% offset
h = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
histogram(o(:,1));
xlabel('offset dim 1 (um)');
x = xlim;
del = x(2)-x(1);

subplot(2,3,2);
histogram(o(:,2));
xlabel('offset dim 2 (um)');
x = xlim;
ndel = x(2)-x(1);
if ndel > del
    del = ndel;
end

subplot(2,3,3);
histogram(o(:,3));
xlabel('offset dim 3 (um)');
x = xlim;
ndel = x(2)-x(1);
if ndel > del
    del = ndel;
end

%
%
%

subplot(2,3,1);
x = xlim;
ndel = x(2)-x(1);
m = (del-ndel)/2;
xlim([x(1)-m x(2)+m]);

subplot(2,3,2);
x = xlim;
ndel = x(2)-x(1);
m = (del-ndel)/2;
xlim([x(1)-m x(2)+m]);

subplot(2,3,3);
x = xlim;
ndel = x(2)-x(1);
m = (del-ndel)/2;
xlim([x(1)-m x(2)+m]);

%
%
%
subplot(2,3,4);
histogram(r(:,1));
xlabel('rotation dim 1 (radians)');
x = xlim;
del = x(2)-x(1);

subplot(2,3,5);
histogram(r(:,2));
xlabel('rotation dim 2 (radians)');
x = xlim;
ndel = x(2)-x(1);
if ndel > del
    del = ndel;
end

subplot(2,3,6);
histogram(r(:,3));
xlabel('rotation dim 3 (radians)');
x = xlim;
ndel = x(2)-x(1);
if ndel > del
    del = ndel;
end


%
%
%

subplot(2,3,4);
x = xlim;
ndel = x(2)-x(1);
m = (del-ndel)/2;
xlim([x(1)-m x(2)+m]);

subplot(2,3,5);
x = xlim;
ndel = x(2)-x(1);
m = (del-ndel)/2;
xlim([x(1)-m x(2)+m]);

subplot(2,3,6);
x = xlim;
ndel = x(2)-x(1);
m = (del-ndel)/2;
xlim([x(1)-m x(2)+m]);

%
%
%

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Repeatability of trajectory of simulated annealing (N = %d)',a(1));
text(0.25,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/histograms.png'];
print(h,str,'-dpng');


end






%%
function plot_raw (oI, oF, rI, rF, out, ppath)

a = size(oI);
% offset
h = figure;

subplot(1,3,1);
for i=1:a(1)
plot( [0 1], [oI(i,1) oF(i,1)], 'k-o' );
hold on;
end
hold off;
ylabel('offset in dim 1 (um)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);

subplot(1,3,2);
for i=1:a(1)
plot( [0 1], [oI(i,2) oF(i,2)], 'k-o' );
hold on;
end
hold off;
ylabel('offset in dim 2 (um)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
for i=1:a(1)
plot( [0 1], [oI(i,3) oF(i,3)], 'k-o' );
hold on;
end
hold off;
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

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Repeatability of trajectory of simulated annealing (N = %d)',a(1));
text(0.25,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_offset.png'];
print(h,str,'-dpng');


% rotation
h = figure;

subplot(1,3,1);
for i=1:a(1)
plot( [0 1], [rI(i,1) rF(i,1)], 'k-o' );
hold on;
end
hold off;
ylabel('rotation around dim 1 (radians)');
set(gca,'xtick',[]);
y = ylim;
del = y(2)-y(1);


subplot(1,3,2);
for i=1:a(1)
plot( [0 1], [rI(i,2) rF(i,2)], 'k-o' );
hold on;
end
hold off;
ylabel('rotation around dim 2 (radians)');
set(gca,'xtick',[]);
y = ylim;
ndel = y(2)-y(1);
if ndel > del
    del = ndel;
end

subplot(1,3,3);
for i=1:a(1)
plot( [0 1], [rI(i,3) rF(i,3)], 'k-o' );
hold on;
end
hold off;
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

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = sprintf('Repeatability of trajectory of simulated annealing (N = %d)',a(1));
text(0.25,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

str = [ppath '/distribution_rotation.png'];
print(h,str,'-dpng');

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
end

function out = get_match (str, content)
expr = ['[^\n]* ' str ':[^\n]*'];
match = regexp(content,expr,'match');
s = regexp(match,'\[.*\]','match');
out = s{1}{1};
end
