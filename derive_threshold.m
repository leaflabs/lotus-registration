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

% if i == 1
%   % i equaling 1 leads to a bug in later calls
%   % just use the next largest value in the cdf
%   % fixes a problem with large sparse images where < 10^-4
%   % fraction of the pixels are non-zero
%   i = 2
% end

end
