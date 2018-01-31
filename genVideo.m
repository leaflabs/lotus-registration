%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%clear all; 
%close all; 
%clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%% Load Source Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data3DPath = '/om/user/jkinney/DLFM/20171124/registration/movie/';
data3DFiles = dir([data3DPath filePattern])
%data3DFiles = dir([data3DPath '*FrameNumber*.mat'])

f =[data3DPath data3DFiles(2).name]
XguessSAVE1 = loadFile(f,eghbal);
% load(f);
% whos
% if eghbal
%     XguessSAVE1 = im2uint8(Xvolume);
%     whos
%     %XguessSAVE1 = Xvolume;
%     clear Xvolume;
% end
movieSize = size(XguessSAVE1);
movieSize(4) = length(data3DFiles);
disp(['Data size :  ' num2str(movieSize(1)) ' x ' num2str(movieSize(2)) ' x ' num2str(movieSize(3)) ' x ' num2str(movieSize(4)) ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nnum = 15;

useInterpolation = true;
GPUinterp = false;
zExpandRatio = 8;

videoRange = 2:length(data3DFiles);

gapMIP = 50;
showVideo = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
singleFrame = XguessSAVE1;
volumeResolution = size(singleFrame);
volumeResolutionRender = size(singleFrame);
%h1 = figure; imshowMIPorth(singleFrame, zExpandRatio, 200, 20, sixteen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure('Visible','on');
iptsetpref('ImshowBorder','tight');

k = 1;
for t=videoRange
    tic;
    %load([data3DPath data3DFiles(t).name]);
    f = [data3DPath data3DFiles(t).name]
    XguessSAVE1 = loadFile(f,eghbal);
    singleFrame = XguessSAVE1;
    imshowMIPorth(singleFrame, zExpandRatio, 200, 20);  
    VideoMIP(k) = getframe(h1);
    ttime = toc;
    disp(['Redering frame ' num2str(k) ' | ' num2str(length(videoRange)) ', took ' num2str(ttime) ' secs']);
    k = k + 1;
end

%timestamp = datestr(datetime('now'));
timestamp = datestr(datetime('now'),'yyyymmdd_HHMMSS');
outputFileName = [dataOutFileName '_' timestamp];
%v = VideoWriter([data3DPath outputFileName],'Uncompressed AVI');
v = VideoWriter([savePath outputFileName],'Uncompressed AVI');
open(v);
writeVideo(v,VideoMIP);
close(v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = loadFile (f, eghbal)
load(f);
if eghbal
    out = Xvolume;
else
    out = XguessSAVE1;
end
end
