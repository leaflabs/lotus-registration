function save_vol (A, savePath, outFile)
imwrite( squeeze(A(:,:,1)), [savePath outFile]);
for k = 2:size(A,3)
    imwrite(squeeze(A(:,:,k)),  [savePath outFile], 'WriteMode', 'append');
end
end