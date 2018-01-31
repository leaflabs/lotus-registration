% from /om/user/jkinney/DLFM/worm/20171215/interpolate/Recon3D_solver_1_FrameNumber_0001__20180123_110051.log
%
% trans = [45.54140 30.9580 52.98630],
% rot = [-1.519490       00       00] 
% final offset = [0.3213920 -8.609540 1.986310]
% centroid: [45.2200 39.5675 51]
% clip: [450 250 150 350 0 0]
%
% from /om/user/jkinney/DLFM/worm/20171215/interpolate/Recon3D_solver_1_FrameNumber_1000__20180123_110443.log
%
% trans = [37.56430 48.03390 55.08210], 
% rot = [-1.517870       00       00]  
% final offset = [-7.655650 8.46640 4.082070]
% centroid: [45.2200 39.5675 51]
% clip: [450 250 150 350 0 0]


m = 1;
n = 1000;
trans1 = [45.54140 30.9580 52.98630];
trans2 = [37.56430 48.03390 55.08210];
rot1 = [-1.519490       00       00];
rot2 = [-1.517870       00       00];

for i=[1 1000]
    clearvars -except i m n trans1 trans2 rot1 rot2;
    close all;
    %rpath = '/om/user/jkinney/DLFM';
    rpath = '/Users/justin/Desktop/DLFM';
    addpath(rpath);
%     dpath = '/om/project/boyden/DualLensLightField/12_15_17_40X_worm/video 1';
%     inputFilePath1 = [dpath '/horizontal/Reconstructed/'];
%     inputFilePath2 = [dpath '/vertical/Reconstructed/'];
%     savePath       = [rpath '/worm/20171215/registration/'];
inputFilePath1 = [rpath '/worm/20171215/horizontal/Reconstructed/'];
inputFilePath2 = [rpath '/worm/20171215/vertical/Reconstructed/'];
savePath       = [rpath '/worm/20171215/registration/'];

    clip    = [450 250 150 350 0 0];
    voxel_x = 0.323/2; % um
    voxel_y = 0.323/2; % um
    voxel_z = 2.0; % um

    rapid = true;
    lfdisplay = false;
    centroid = [45.2200 39.5675 51] ...
        + [clip(1)*voxel_y clip(3)*voxel_x 0];
    rot = rot1 + (rot1-rot2) * (i-m)/(n-m);
    trans = trans1 + (trans1-trans2) * (i-m)/(n-m);
    trans = trans + [clip(1)*voxel_y clip(3)*voxel_x 0];
    inputFileName = {sprintf('Recon3D_solver_1_FrameNumber_%04d.mat',i)};
    register
end
