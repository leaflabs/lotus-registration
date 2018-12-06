function psfVals = plot_psf_cropped (config, handle, csec, dim, str, Buf)
subplot(handle);
a = max(csec,[],dim);
psfVals = single(a)/max(single(a(Buf:length(a)-Buf+1)));
vec = [1:length(psfVals)]*config.pixel;
plot(vec,psfVals,'*-');
%set(gca, 'YDir', 'reverse');
%thisdim = get(gca,'Position');
%set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
%set(gca,'XTickLabel',{});
%b = size(csec);
%xlim([1 b(2)]);
ylabel('Normalized intensity');
ylim([0 1]);
xlabel(str);
end