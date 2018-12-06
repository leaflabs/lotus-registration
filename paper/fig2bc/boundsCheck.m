function out = boundsCheck (cond, val, bound, mystr, a)
if cond
    msg = sprintf('Warning. %s pixel %d exceeds bound %d and is outside slice (%d,%d).', ...
        mystr,val, bound, a(1),a(2));
    disp(msg);
    %out = bound;
    %msg = sprintf('%s pixel set to %d.', mystr, out);
    msg = sprintf('Skipping bead.');
    disp(msg);
    out = -1;
else
    out = val;
end
end