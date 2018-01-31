function mat2tif( f, o )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load(f);

if ~exist('XguessSAVE1','var')
    keyboard
end

save_vol(XguessSAVE1, '', o);

end


function save_vol (A, savePath, outFile)
imwrite( squeeze(A(:,:,1)), [savePath outFile]);
for k = 2:size(A,3),
    imwrite(squeeze(A(:,:,k)),  [savePath outFile], 'WriteMode', 'append');
end
end