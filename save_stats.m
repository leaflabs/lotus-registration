function save_stats (param)

prefix = sprintf('%s%s_',param.savePath,param.timestamp);


% plot Pvec
h = figure;
plot(1:numel(param.Pvec),param.Pvec);
xlabel('iteration');
ylabel('probability of allowing a decrease in MI');
title('simulated annealing schedule');
if 0>1
    fname = sprintf('%s_p.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_p.png',prefix);
    save_plot(h, str);
end


% plot Tvec
h = figure;
plot(1:numel(param.Tvec),param.Tvec);
xlabel('iteration');
ylabel('Temperature (-)');
title('simulated annealing schedule');
if 0>1
    fname = sprintf('%s_p.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_T.png',prefix);
    save_plot(h, str);
end

% plot MIvec
h = figure;
plot(1:numel(param.MIvec),log10(param.MIvec));
hold on;
plot(xlim,log10(param.bestMI)*ones(2,1),'--r');
hold off;
xlabel('iteration');
ylabel('log (mutual information)');
title('evolution of mutual information during simulated annealing');
if 0>1
    fname = sprintf('%s_MI.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_MI.png',prefix);
    save_plot(h, str);
end

% plot CDFvec
f = figure;
plot(1:numel(param.cdfvec),param.cdfvec);
hold on;
plot(xlim,[1 1],'--r');
hold off;
xlabel('iteration');
ylabel('fraction of null distribution');
title('fraction of null distribution less than current mutual informaton');
if 0>1
    fname = sprintf('%s_frac_null.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_frac_null.png',prefix);
    save_plot(f, str);
end

% plot amplitude of 3D random walk
h = figure;
a = size(param.Dvec);
subplot(1,3,1);
plot(1:a(1),param.Dvec(:,1));
xlabel('iteration');
ylabel('random walk in dim 1 (um)');
hold on;
plot(1,param.Dvec(1,1),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,2);
plot(1:a(1),param.Dvec(:,2));
xlabel('iteration');
ylabel('random walk in dim 2 (um)');
hold on;
plot(1,param.Dvec(1,2),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,3);
plot(1:a(1),param.Dvec(:,3));
xlabel('iteration');
ylabel('random walk in dim 3 (um)');
hold on;
plot(1,param.Dvec(1,3),'ro');
hold off;
xlim([1 a(1)]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'Attempted translations';
text(0.4,0.97,str,'FontSize',12,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_d.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_d.png',prefix);
    save_plot(h, str);
end

% plot offset
h = figure;
a = size(param.offsetvec);
subplot(1,3,1);
plot(1:a(1),param.offsetvec(:,1));
xlabel('iteration');
ylabel('offset in dim 1 (um)');
hold on;
plot(1,param.offsetvec(1,1),'ro');
plot(a(1),param.bestd(1),'gs');
hold off;
xlim([1 a(1)]);

subplot(1,3,2);
plot(1:a(1),param.offsetvec(:,2));
xlabel('iteration');
ylabel('offset in dim 2 (um)');
hold on;
plot(1,param.offsetvec(1,2),'ro');
plot(a(1),param.bestd(2),'gs');
hold off;
xlim([1 a(1)]);

subplot(1,3,3);
plot(1:a(1),param.offsetvec(:,3));
xlabel('iteration');
ylabel('offset in dim 3 (um)');
hold on;
plot(1,param.offsetvec(1,3),'ro');
plot(a(1),param.bestd(3),'gs');
hold off;
xlim([1 a(1)]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'trajectory of simulated annealing (Green squares are best of null)';
text(0.4,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_offset.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_offset.png',prefix);
    save_plot(h, str);
end

% plot rot
h = figure;
plot_rot(h,1,param);
plot_rot(h,2,param);
plot_rot(h,3,param);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
%str = 'trajectory of simulated annealing (RHS axis is degrees, green squares are best of null.)';
str = 'trajectory of simulated annealing (RHS axis is degrees)';
text(0.1,0.98,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_rot.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_rot.png',prefix);
    save_plot(h, str);
end

end
