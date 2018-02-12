%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%clear all; 
%close all; 
%clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%% Load Source Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nLooking for LFM1 volume files with this path and pattern:\n%s\n\n',LFM1_Path);
LFM1_Files = dir(LFM1_Path);
fprintf('Found %d files\n\n',length(LFM1_Files));

% get size of reconstructed volumes
% arbitrarily pick the first volume to measure
LFM1_movieSize = get_size ([LFM1_Files(1).folder '/' LFM1_Files(1).name]);

LFM1_movieSize(4) = length(LFM1_Files);
disp(['Data size :  ' num2str(LFM1_movieSize(1)) ' x ' num2str(LFM1_movieSize(2))...
    ' x ' num2str(LFM1_movieSize(3)) ' x ' num2str(LFM1_movieSize(4)) ]);





fprintf('\nLooking for LFM2 volume files with this path and pattern:\n%s\n\n',LFM2_Path);
LFM2_Files = dir(LFM2_Path);
fprintf('Found %d files\n\n',length(LFM2_Files));

% get size of reconstructed volumes
% arbitrarily pick the first volume to measure
LFM2_movieSize = get_size ([LFM2_Files(1).folder '/' LFM2_Files(1).name]);

LFM2_movieSize(4) = length(LFM2_Files);
disp(['Data size :  ' num2str(LFM2_movieSize(1)) ' x ' num2str(LFM2_movieSize(2))...
    ' x ' num2str(LFM2_movieSize(3)) ' x ' num2str(LFM2_movieSize(4)) ]);





fprintf('\nLooking for DLFM volume files with this path and pattern:\n%s\n\n',DLFM_Path);
DLFM_Files = dir(DLFM_Path);
fprintf('Found %d files\n\n',length(DLFM_Files));

% get size of reconstructed volumes
% arbitrarily pick the first volume to measure
DLFM_movieSize = get_size ([DLFM_Files(1).folder '/' DLFM_Files(1).name]);

DLFM_movieSize(4) = length(DLFM_Files);
disp(['Data size :  ' num2str(DLFM_movieSize(1)) ' x ' num2str(DLFM_movieSize(2))...
    ' x ' num2str(DLFM_movieSize(3)) ' x ' num2str(DLFM_movieSize(4)) ]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nnum = 15;

useInterpolation = true;
GPUinterp = false;
% zExpandRatio = 8;
zExpandRatio = round(4 / 0.323);
zExpandRatioDLFM = zExpandRatio / 8;
moviePixelVert = 1000;

videoRange = 1:length(LFM1_Files);

gapMIP = 20;  % size of gap in units of pixels
gapVal = 200; % value of gap pixels
gapMETA = 200;
showVideo = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%volumeResolution = LFM1_movieSize;
%volumeResolutionRender = LFM1_movieSize;
%h1 = figure; imshowMIPorth(singleFrame, zExpandRatio, 200, 20, sixteen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure('Visible','on');
iptsetpref('ImshowBorder','tight');

k = 1;
for t=videoRange
    tic;
    LFM1 = loadFile([LFM1_Files(t).folder '/' LFM1_Files(t).name]);
    LFM2 = loadFile([LFM2_Files(t).folder '/' LFM2_Files(t).name]);
    DLFM = loadFile([DLFM_Files(t).folder '/' DLFM_Files(t).name]);
    imshowMIPorth(LFM1, LFM2, DLFM, zExpandRatio, zExpandRatioDLFM, gapVal, gapMIP, gapMETA, moviePixelVert);
    VideoMIP(k) = getframe(h1);
    ttime = toc;
    disp(['Redering frame ' num2str(k) ' | ' num2str(length(videoRange)) ', took ' num2str(ttime) ' secs']);
    k = k + 1;
end

timestamp = datestr(datetime('now'),'yyyymmdd_HHMMSS');
outputFileName = [dataOutFileName '_' timestamp];
v = VideoWriter([savePath outputFileName],'Uncompressed AVI');
open(v);
writeVideo(v,VideoMIP);
close(v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = loadFile (f)
load(f);
if exist('Xvolume') == 1
    out = Xvolume;
elseif exist('XguessSAVE1') == 1
    out = XguessSAVE1;
else
    disp('WTF?');
    whos
    keyboard
end
end

function out = get_size (f)
load(f);
if exist('Xvolume') == 1
    out = size(Xvolume);
elseif exist('XguessSAVE1') == 1
    out = size(XguessSAVE1);
else
    disp('WTF?');
    whos
    keyboard
end
end



function IMGout = imshowMIPorth(LFM1, LFM2, DLFM, zExpansionFactor, zExpansionFactorDLFM, gapVal, gapMIP, gapMETA, moviePixelVert)
[xsize1, ysize1] = get_MIPsize (LFM1, zExpansionFactor, gapMIP);
[xsize2, ysize2] = get_MIPsize (LFM2, zExpansionFactor, gapMIP);
[xsize3, ysize3] = get_MIPsize (DLFM, zExpansionFactorDLFM, gapMIP);
x = (xsize1 + gapMETA) + (xsize2 + gapMETA) + xsize3;
%y = max( [ ysize1 ysize2 ysize3] );

IMGout = gapVal*ones(moviePixelVert, x, 'uint16');

origin = 1;
s = size(LFM1);
IMGout = metaMIPs(LFM1, s(2), s(1), zExpansionFactor*s(3), gapMIP, origin, IMGout);
origin = (xsize1 + gapMETA) + 1;
s = size(LFM2);
IMGout = metaMIPs(LFM2, s(2), s(1), zExpansionFactor*s(3), gapMIP, origin, IMGout);
origin = (xsize1 + gapMETA) + (xsize2 + gapMETA) + 1;
s = size(DLFM);
IMGout = metaMIPs(DLFM, s(2), s(1), zExpansionFactorDLFM*s(3), gapMIP, origin, IMGout);

imshow(IMGout);

end


function [xsize,ysize] = get_MIPsize (LFM, zExpansionFactor, gapMIP)
s = size(LFM);
xsize = s(2) + zExpansionFactor*s(3) + gapMIP;
ysize = s(1) + zExpansionFactor*s(3) + gapMIP;
end



function IMGout = metaMIPs(LFM, sizex, sizey, sizez, gapMIP, origin, IMGout)

if ~isa(LFM,'uint16')
    disp('Warning input data is not 16 bit.');
    keyboard
end
IMGin1 =           squeeze(max(LFM,[],3));
IMGin2 = imresize( squeeze(max(LFM,[],1)), [sizex sizez ]);
IMGin3 = imresize( squeeze(max(LFM,[],2)), [sizey sizez ]);

xsize = sizex + sizez + gapMIP;
ysize = sizey + sizez + gapMIP;

IMGout( (1:sizey), ((origin-1)+1:(origin-1)+sizex) ) = IMGin1;
IMGout( (1:sizey), ((origin-1)+sizex+gapMIP+1:(origin-1)+xsize) ) = IMGin3;
IMGout( (sizey+gapMIP+1:ysize), ((origin-1)+1:(origin-1)+sizex) ) = IMGin2';

end



% function IMGout = combineMIPs(IMGin1, IMGin2, IMGin3, sizex, sizey, sizez, gapVal, gapMIP)
% 
% xsize = sizex + sizez + gapMIP;
% ysize = sizey + sizez + gapMIP;
% 
% IMGout = gapVal*ones(ysize, xsize, 'uint16');
% 
% IMGout( (1:sizey), (1:sizex) ) = IMGin1;
% IMGout( (1:sizey), (sizex+gapMIP+1:xsize) ) = IMGin3;
% IMGout( (sizey+gapMIP+1:ysize), (1:sizex) ) = IMGin2';
% end
