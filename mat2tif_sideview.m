function mat2tif_sideview( f, resize )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
o = [f(1:end-4) '_sideview.tif'];
load(f);
if exist('XguessSAVE1','var')
    save_vol(XguessSAVE1, '', o, resize);
elseif exist('Xvolume','var')
    save_vol(Xvolume, '', o, resize);
else
    keyboard
end
end


function save_vol (A, savePath, outFile, resize)
if isempty(resize)
    imwrite( squeeze(A(:,1,:)), [savePath outFile]);
else
    B = squeeze(A(:,1,:));
    a = size(B);
    imwrite( imresize( B, [a(1) 8*a(2)] ), [savePath outFile]);
end

for k = 2:size(A,2)
    if isempty(resize)
        imwrite(squeeze(A(:,k,:)),  [savePath outFile], 'WriteMode', 'append');
    else
        B = squeeze(A(:,k,:));
        a = size(B);
        imwrite( imresize( B, [a(1) 8*a(2)] ),  [savePath outFile], 'WriteMode', 'append');
    end
end
end