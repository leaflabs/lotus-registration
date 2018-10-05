function mat2tif( f )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
o = [f(1:end-4) '.tif'];
load(f);
if exist('XguessSAVE1','var')
    save_vol(XguessSAVE1, '', o);
elseif exist('Xvolume','var')
    save_vol(Xvolume, '', o);
else
    keyboard
end
end


function save_vol (A, savePath, outFile)
imwrite( squeeze(A(:,:,1)), [savePath outFile]);
for k = 2:size(A,3),
    imwrite(squeeze(A(:,:,k)),  [savePath outFile], 'WriteMode', 'append');
end
end