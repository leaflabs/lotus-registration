function i = derive_threshold (cdf, param)
% look for large jump in cdf
dif = cdf(2:end)-cdf(1:end-1);
sdif = flipud(sortrows([dif' [1:length(dif)]'],1));
if sdif(1,1) > param.dynamic_range_thresh * sdif(2,1)
    i = sdif(1,2)+1;
else
    % if not found
    i = find(cdf>param.pop_thresh,1); % keep only first instance
end
end