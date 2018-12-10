function print_fwhm_ratio( bs, ws, str )
tmp1 = mean(bs(ws,6) ./ bs(ws,8));
tmp2 = mean(bs(ws,7) ./ bs(ws,8));
fprintf('%s mean [1:3__%2.1f, 2:3__%2.1f]\n',str,tmp1,tmp2);
tmp1 = std(bs(ws,6) ./ bs(ws,8));
tmp2 = std(bs(ws,7) ./ bs(ws,8));
fprintf('     stddev [1:3__%2.1f, 2:3__%2.1f]\n',tmp1,tmp2);
end

