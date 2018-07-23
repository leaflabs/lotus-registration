function save_2d_contour_plots (LFM1, LFM2, pos1, pos2, new, param, str)
colors = {[1 0 0]
    [0.8 0.8 0.8]
    [0 0 1]
    [0.4 0.4 0.4]};

% projections for new
% convert new from microns to pixels
a = size(new);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_new = ceil(new.*scale);
ns = max(abc_new);

% create 2d projections of new in pixel space
onetwo_new = zeros(ns(1),ns(2));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,2) > 0
        onetwo_new(abc_new(i,1),abc_new(i,2)) ...
            = max([onetwo_new(abc_new(i,1),abc_new(i,2)) ...
            double(LFM2(param.index2(i)))]);
    end
end
onethree_new = zeros(ns(1),ns(3));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,3) > 0
        onethree_new(abc_new(i,1),abc_new(i,3)) ...
            = max([onethree_new(abc_new(i,1),abc_new(i,3))...
            double(LFM2(param.index2(i)))]);
    end
end
threetwo_new = zeros(ns(3),ns(2));
for i=1:length(new)
    if abc_new(i,3) > 0 && abc_new(i,2) > 0
        threetwo_new(abc_new(i,3),abc_new(i,2))...
            = max([threetwo_new(abc_new(i,3),abc_new(i,2))...
            double(LFM2(param.index2(i)))]);
    end
end


% projections for pos1
% convert new from microns to pixels
a = size(pos1);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_pos1 = ceil(pos1.*scale);
ns = max(abc_pos1);

% create 2d projections of new in pixel space
onetwo_pos1 = zeros(ns(1),ns(2));
for i=1:length(pos1)
    if abc_pos1(i,1) > 0 && abc_pos1(i,2) > 0
        onetwo_pos1(abc_pos1(i,1),abc_pos1(i,2))...
            = max([onetwo_pos1(abc_pos1(i,1),abc_pos1(i,2))...
            double(LFM1(param.index1(i)))]);
    end
end
onethree_pos1 = zeros(ns(1),ns(3));
for i=1:length(pos1)
    if abc_pos1(i,1) > 0 && abc_pos1(i,3) > 0
        onethree_pos1(abc_pos1(i,1),abc_pos1(i,3))...
            = max([onethree_pos1(abc_pos1(i,1),abc_pos1(i,3))...
            double(LFM1(param.index1(i)))]);
    end
end
threetwo_pos1 = zeros(ns(3),ns(2));
for i=1:length(pos1)
    if abc_pos1(i,3) > 0 && abc_pos1(i,2) > 0
        threetwo_pos1(abc_pos1(i,3),abc_pos1(i,2))...
            = max([threetwo_pos1(abc_pos1(i,3),abc_pos1(i,2))...
            double(LFM1(param.index1(i)))]);
    end
end

% projections for pos2
% convert new from microns to pixels
a = size(pos2);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_pos2 = ceil(pos2.*scale);
ns = max(abc_pos2);

% create 2d projections of new in pixel space
onetwo_pos2 = zeros(ns(1),ns(2));
for i=1:length(pos2)
    if abc_pos2(i,1) > 0 && abc_pos2(i,2) > 0
        onetwo_pos2(abc_pos2(i,1),abc_pos2(i,2)) ...
            = max([onetwo_pos2(abc_pos2(i,1),abc_pos2(i,2)) ...
            double(LFM2(param.index2(i)))]);
    end
end
onethree_pos2 = zeros(ns(1),ns(3));
for i=1:length(pos2)
    if abc_pos2(i,2) > 0 && abc_pos2(i,3) > 0
        onethree_pos2(abc_pos2(i,1),abc_pos2(i,3)) ...
            = max([onethree_pos2(abc_pos2(i,1),abc_pos2(i,3)) ...
            double(LFM2(param.index2(i)))]);
    end
end
threetwo_pos2 = zeros(ns(3),ns(2));
for i=1:length(pos2)
    if abc_pos2(i,3) > 0 && abc_pos2(i,2) > 0
        threetwo_pos2(abc_pos2(i,3),abc_pos2(i,2)) ...
            = max([threetwo_pos2(abc_pos2(i,3),abc_pos2(i,2)) ...
            double(LFM2(param.index2(i)))]);
    end
end

% centroid in pixels
c = [param.centroid(1)/param.voxel_y...
    param.centroid(2)/param.voxel_x...
    param.centroid(3)/param.voxel_z];

f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
hold on;
a = size(onethree_new);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onethree_new,v,'Color',colors{1});
a = size(onethree_pos1);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,onethree_pos1,v,'Color',colors{3});
a = size(onethree_pos2);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onethree_pos2,v,'Color',colors{2});
xlabel('dim three [pixels]');
ylabel('dim one [pixels]');
plot(c(3),c(1),'g*')
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');
set(gca,'XDir','reverse');


subplot(2,2,2);
hold on;
a = size(onetwo_new);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onetwo_new,v,'Color',colors{1});
a = size(onetwo_pos1);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,onetwo_pos1,v,'Color',colors{3});
a = size(onetwo_pos2);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onetwo_pos2,v,'Color',colors{2});
xlabel('dim two [pixels]');
ylabel('dim one [pixels]');
plot(c(2),c(1),'g*')
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

subplot(2,2,4);
hold on;
a = size(threetwo_new);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,threetwo_new,v,'Color',colors{1});
a = size(threetwo_pos1);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,threetwo_pos1,v,'Color',colors{3});
a = size(threetwo_pos2);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,threetwo_pos2,v,'Color',colors{2});
xlabel('dim two [pixels]');
ylabel('dim three [pixels]');
plot(c(2),c(3),'g*');
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);

text(0.2,0.45,str,'FontSize',12,'Color',[0 0 0] ,'Interpreter','none');

msg = sprintf('LFM1: contour at intensity level %4.0f',param.contour_int1);
text(0.2,0.4,msg,'FontSize',12,'Color',colors{3},'Interpreter','none');

msg = sprintf('LFM2: contour at intensity level %4.0f',param.contour_int2);
text(0.2,0.35,msg,'FontSize',12,'Color',colors{2},'Interpreter','none');

msg = sprintf('LFM2 coarse reg: contour at intensity level %4.0f',param.contour_int2);
text(0.2,0.3,msg,'FontSize',12,'Color',colors{1},'Interpreter','none');

text(0.2,0.25,'centroid of LFM2 coarse reg','FontSize',12,'Color',[0 1 0],'Interpreter','none');

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
save_plot(f, str);
end