function IMGout = combineMIPs(IMGin1, IMGin2, IMGin3, sizex, sizey, sizez, gapVal, gapMIP)

xsize = sizex + sizez + gapMIP;
ysize = sizey + sizez + gapMIP;

IMGout = gapVal*ones(ysize, xsize, 'uint16');

IMGout( (1:sizey), (1:sizex) ) = IMGin1;
IMGout( (1:sizey), (sizex+gapMIP+1:xsize) ) = IMGin3;
IMGout( (sizey+gapMIP+1:ysize), (1:sizex) ) = IMGin2';

