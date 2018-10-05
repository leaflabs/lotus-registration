function xyz = getXYZbead (vol, outpath, prefix)
xy = squeeze(sum(vol,3));
xyz = [];

f = figure(1);
dynamic_range = round(log2(max(max(xy))));
imagesc(xy,[0 2^dynamic_range]);
daspect([1,1,1]);
xlabel('x');
ylabel('y');
colorbar();
xycoords = [];
hold on;
while (1)
    a = round(ginput(1));
    if a(1)<50 & a(2) <50
        break;
    end
    xycoords = [xycoords;a];
    plot(a(1),a(2),'k*');
    row = a(2);
    col = a(1);
    z = get_z( row, col, vol );
    xyz = [xyz;[a z]];
    msg = sprintf('bead found at row = %d, col = %d, z = %d',row,col,z);
    disp(msg);
end
hold off;

fout = [outpath prefix '_bead_xyz.mat'];
save(fout,'xyz');
fprintf('\n# Output file = %s\n',fout);

fout = [outpath prefix '_bead_xy.png']
print(f,fout,'-dpng');
fprintf('\n# Output file = %s\n',fout);
end


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