function fwhm = plot_psf_cropped (handle, csec, dim, str, scale)
a = max(csec,[],dim);
psfVals = double(a)/max(double(a));
vec = [1:length(psfVals)]*scale;
[y,i] = max(a);
vec_c = vec - vec(i);
a = size(psfVals);
if a(2)==1
    psfVals = psfVals';
end
[f, gof] = fit(vec_c',psfVals','gauss1');
%
subplot(handle);
p = plot(f,vec_c,psfVals);
p(1).MarkerSize = 50;
p(2).LineWidth = 8;
ylabel('Intensity (a.u.)');
ylim([0 1]);
%
xlabel(str);
ax=gca;
legend(ax,'off');
ax.YTick = [0 1];
ax.FontSize = 40;
%
hwhm = f.b1 + f.c1 * sqrt(-log(0.5/f.a1));
fwhm = 2 * hwhm;
%
L = 10;
H = 30;
if fwhm < 10
    V = L;
    ax.XTick = [-10 -5 0 5 10];
else
    V = H;
    ax.XTick = [-30 -15 0 15 30];
end
xlim([-V V]);
%
str = sprintf('fwhm = %2.1f um',fwhm);
text(-V,0.5,str,'FontSize',40,'Color',[0 0 0],'Interpreter','none');
end