function csec_cropped = cropImag(csec, crop_index)
csec_cropped = csec(crop_index(1):crop_index(2),crop_index(3):crop_index(4));
a = size(csec);
msg = sprintf('Orig slice is %d by %d',a(1),a(2));
disp(msg);
b = size(csec_cropped);
msg = sprintf('Cropped slice is %d by %d',b(1),b(2));
disp(msg);
end