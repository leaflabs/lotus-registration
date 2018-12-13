function plot_fwhm ( config, bead_struct )
ms = 20;
c = [0.5 0.5 0.5];

wl_dlfm = config.whitelist_dlfm;
wl_lfm1 = config.whitelist_lfm1;
wl_lfm2 = config.whitelist_lfm2;

max_dim1 = max ([max(bead_struct.dlfm(wl_dlfm,6)) ...
    max(bead_struct.lfm1(wl_lfm1,6)) ...
    max(bead_struct.lfm2(wl_lfm2,6))] );

max_dim2 = max ([max(bead_struct.dlfm(wl_dlfm,7)) ...
    max(bead_struct.lfm1(wl_lfm1,7)) ...
    max(bead_struct.lfm2(wl_lfm2,7))] );

max_dim3 = max ([max(bead_struct.dlfm(wl_dlfm,8)) ...
    max(bead_struct.lfm1(wl_lfm1,8)) ...
    max(bead_struct.lfm2(wl_lfm2,8))] );

mystr = {'DLFM','LFM1','LFM2'};

f = figure();
set(f,'Position',config.figpos);

subplot(1,2,1);
hold on;
plot(bead_struct.dlfm(wl_dlfm,8),bead_struct.dlfm(wl_dlfm,6),...
    '.','MarkerSize',ms);
plot(bead_struct.lfm1(wl_lfm1,8),bead_struct.lfm1(wl_lfm1,6),...
    '<','MarkerSize',ms,'Color',c);
plot(bead_struct.lfm2(wl_lfm2,8),bead_struct.lfm2(wl_lfm2,6),...
    '>','MarkerSize',ms,'Color',c);
p = plot([0 40],[0 40],'k--');
hold off;
legend(mystr);
axis equal;
xlabel('fwhm in dim 3 or z (um)');
ylabel('fwhm in dim 1 or y (um)');
%
xlim([0 max_dim3+10]);
ylim([0 max_dim1+10]);
ax=gca;
ax.FontSize = 40;

subplot(1,2,2);
hold on;
plot(bead_struct.dlfm(wl_dlfm,7),bead_struct.dlfm(wl_dlfm,6),...
'.','MarkerSize',ms);
plot(bead_struct.lfm1(wl_lfm1,7),bead_struct.lfm1(wl_lfm1,6),...
'<','MarkerSize',ms,'Color',c);
plot(bead_struct.lfm2(wl_lfm2,7),bead_struct.lfm2(wl_lfm2,6),...
'>','MarkerSize',ms,'Color',c);
plot([0 40],[0 40],'k--');
axis equal;
xlabel('fwhm in dim 2 or x (um)');
ylabel('fwhm in dim 1 or y (um)');
hold off;
legend(mystr);
v = max(max_dim1,max_dim2);
xlim([0 v+10]);
ylim([0 v+10]);
ax=gca;
ax.FontSize = 40;

str=sprintf('%s%s%s_fwhm_summary.png',config.outpath,config.label,config.fname(1:end-4));
print(f,str,'-dpng');
str=sprintf('%s%s%s_fwhm_summary.eps',config.outpath,config.label,config.fname(1:end-4));
print(f,str,'-depsc');

end

