function print_ms( bs, wl, str )
tmp = mean(bs(wl,6:8));
fprintf('%s mean [%2.1f %2.1f %2.1f]\n',str, tmp(1),tmp(2),tmp(3));
tmp = std(bs(wl,6:8));
fprintf('     stddev [%2.1f %2.1f %2.1f]\n',tmp(1),tmp(2),tmp(3));

end

