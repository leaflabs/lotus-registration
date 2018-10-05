function out = print_param (param)
a=param.trans;
b=param.rot;
out = sprintf('trans = [%6.6f %6.6f %6.6f], rot = [%7.6f %7.6f %7.6f]',a(1),a(2),a(3),b(1),b(2),b(3));
end
