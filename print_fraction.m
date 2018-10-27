function print_fraction (index, LFM, str)
i = 0;
for j=1:length(index)
    i = i + single(LFM(index(j)));
end
iT = single(sum(sum(sum(LFM))));

fi = i/iT;
fv = single(length(index))/single(numel(LFM));
fprintf('%s: %d of %d voxels (%7.3g), %d of %d total intensity (%7.3g)\n',...
    str,length(index),numel(LFM),fv,...
    i,iT,fi);
end
