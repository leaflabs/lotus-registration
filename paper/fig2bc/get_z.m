function z = get_z ( row, col, mymatrix )
a = size(mymatrix);
if row > a(1)
    msg = sprintf('Skipping bead, row %d is outside of image (%d)',row,a(1));
    disp(msg);
    z=-1;
    return;
elseif col > a(2)
    msg = sprintf('Skipping bead, col %d is outside of image (%d)',col,a(2));
    disp(msg);
    z=-1;
    return;
end
mymax = -1;
z = -1;
for i=1:a(3)
    if row<1 | col<1 | i<1
        disp('WTF?');
        keyboard;
    end
    val = mymatrix(row,col,i);
    if val > mymax
        mymax = val;
        z = i;
    end
end
end