function p = calc_boltzman_p (T, param)
val = -param.scale/T;
if val > 0.0
    p = 1.0;
else
    p = exp(val);
end
if p<0.0 || p>1.0
    disp('WTF?');
    keyboard
end
end