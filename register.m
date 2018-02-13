function register (param)

%%
param.timestamp = [param.inputFileName{1}(1:end-4) '__' datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%s%s.log',param.savePath, param.timestamp);
if ~exist(param.savePath,'dir')
    status = mkdir(param.savePath);
    if status == 1
        disp(['Created folder: ' param.savePath]);
    else
        disp(['Error attempting to create folder:' param.savePath]);
        status
        exit;
    end
end
diary(fname)
tic

%% load volumes
f = [param.inputFilePath1 param.inputFileName{1}];
fprintf('\nLoading volume %s\n\n',f);
LFM1 = loadData(f, param);
f = [param.inputFilePath2 param.inputFileName{1}];
fprintf('\nLoading volume %s\n\n',f);
LFM2 = loadData(f, param);
param.voxel_z = param.voxel_z / param.interp;

%% if registration parameters are already known
%  then combine registered volumes
if param.rapid
    % must specify: param.centroid, param.trans, param.rot
    combineVols_iter_dim3 (LFM1, LFM2, param);
    return;
end


%% save TIF versions of input volumes
if param.lfdisplay
    f = param.inputFileName{1};
    outFile = sprintf('%s.tif',f(1:end-4));
    save_vol( LFM1, param.inputFilePath1, outFile);
    save_vol( LFM2, param.inputFilePath2, outFile);
end


%% clip volumes for registration
a = size(LFM1);
if ~isempty(find(param.clip>0))
    fprintf('\nClipping pixels from periphery:\n');
    fprintf('%d \n',param.clip);
    fprintf('\n\n');
    LFM1 = LFM1( 1+param.clip(1):a(1)-param.clip(2),...
        1+param.clip(3):a(2)-param.clip(4),...
        1+param.clip(5):a(3)-param.clip(6));
end
a = size(LFM2);
if ~isempty(find(param.clip>0))
    LFM2 = LFM2( 1+param.clip(1):a(1)-param.clip(2),...
        1+param.clip(3):a(2)-param.clip(4),...
        1+param.clip(5):a(3)-param.clip(6));
end

%% calculate thresholds
param = calculate_thresholds (LFM1, LFM2, param,'cdf_voxel_intensities');

%% voxel positions
param.index1 = find(LFM1>param.threshold1);
print_fraction(param.index1,LFM1,'LFM1');
pos1 = init_pos(param.index1,LFM1,param);
param.index2 = find(LFM2>param.threshold2);
print_fraction(param.index2,LFM2,'LFM2');
%param.size = size(LFM2);
pos2 = init_pos(param.index2,LFM2,param);

%% coarse alignment of LFM2 to LFM1
param.centroid = calc_centroid(LFM2,param);
canonical = translate (pos2, -param.centroid);
rotated = rotate (canonical,param.rot+param.angle);
param.rot = param.rot + param.angle;
%new = translate (rotated, param.centroid+param.offset);
new = translate (rotated, param.centroid);
%param.trans = param.centroid + param.offset;
param.trans = param.centroid;

%% estimate offsets
offsets = estimate_offsets(LFM1, LFM2, new, param);
for i=1:3
    if ~isempty(param.offset{i})
        offsets(i) = param.offset{i};
    end
end
param.offset = offsets;
fprintf('\nestimated offsets = [%f %f %f]\n',offsets(1),offsets(2),offsets(3));
new = translate (new, offsets);
param.trans = param.trans + offsets;

%% plot data
if param.plot
    param = save_1d_max_projections(LFM1, LFM2, pos1, pos2, new, param,...
        '1d_max_projections_presim');
    save_2d_max_projections(LFM1, LFM2, new, param, 0,...
        '2d_max_projections_presim');
    save_2d_max_projections_compact(LFM1, LFM2, new, param, 0,...
        '2d_max_projections_compact_presim');
    save_2d_contour_plots(LFM1, LFM2, pos1, pos2, new, param,...
        'contour_plots_presim');
    drawnow
end

%% null distribution
% for random rotations and offsets (within some range)
% calculate mutual information, MI
[cdf,centers,nullMIvec] = null_distribution (LFM1, LFM2, canonical, param);
param.cdf = cdf;
param.centers = centers;
param.nullMIvec = nullMIvec;

%% print MI of volumes pre- and post-coarse registration
% pre-coarse
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, pos2, LFM2, param);
else
    disp('WTF!');
    keyboard;
end
w = find(param.centers>MI,1); % keep only first instance
if MI > max(nullMIvec)
    fprintf('\npre-coarse MI frac = %.5g\n\n',double(MI)/double(max(nullMIvec)));
elseif ~isempty(w)
    fprintf('\npre-coarse MI frac = %.5g\n\n',param.cdf(w));
else
    disp('WTF?');
    keyboard
end
% post-coarse
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, new, LFM2, param);
else
    disp('WTF!');
    keyboard;
end
w = find(param.centers>MI,1); % keep only first instance
if MI > max(nullMIvec)
    fprintf('\npost-coarse MI frac = %.5g\n\n',double(MI)/double(max(nullMIvec)));
elseif ~isempty(w)
    fprintf('\npost-coarse MI frac = %.5g\n\n',param.cdf(w));
else
    disp('WTF?');
    keyboard
end

%% simulated annealing
[new, param] = simulated_annealing (LFM1, new, LFM2, canonical, param);
if param.plot
    save_stats(param);
end


%% combine registered volumes
if param.savevol
    % would be nice to re-load full unclipped volume here
    % but doing so would jack up the next few plots.
    comb = combineVols (LFM1, LFM2, param);
end

if param.plot
    param = save_1d_max_projections(LFM1, LFM2, pos1, pos2, new, param,...
        '1d_max_projections_postsim');
    save_2d_max_projections(LFM1, LFM2, comb, param, 1, ...
        '2d_max_projections_postsim');
    save_2d_max_projections_compact(LFM1, LFM2, comb, param, 1, ...
        '2d_max_projections_compact_postsim');
    save_2d_contour_plots(LFM1, LFM2, pos1, pos2, new, param, ...
        'contour_plots_postsim');
    drawnow
end

%% save final parameters
%fname = sprintf('%s%s_parameters.mat',param.savePath,param.timestamp);
%save(fname,'param');

param
elapsedTime = toc
diary off;

end


%% functions

function param = calculate_thresholds (LFM1, LFM2, param, str)
% assume 16 bit
edges = linspace(0,2^16,param.N);
centers = (edges(1:end-1)+edges(2:end))/2;
% this works. to check: sum(h_LFM1) == numel(LFM1)
h_LFM1 = histcounts(LFM1,edges);
h_LFM2 = histcounts(LFM2,edges);

% normalize cdf to 1
h_LFM1 = h_LFM1 / sum(sum(sum(h_LFM1)));
h_LFM2 = h_LFM2 / sum(sum(sum(h_LFM2)));

total = 0;
cdf_LFM1 = [];
for i=1:length(h_LFM1)
    total = total + h_LFM1(i);
    cdf_LFM1 = [cdf_LFM1 total];
end
total = 0;
cdf_LFM2 = [];
for i=1:length(h_LFM2)
    total = total + h_LFM2(i);
    cdf_LFM2 = [cdf_LFM2 total];
end

i = derive_threshold (cdf_LFM1, param);
j = derive_threshold (cdf_LFM2, param);

param.threshold1 = centers(i);
param.threshold2 = centers(j);
fprintf('threshold1 = %f\n',param.threshold1);
fprintf('threshold2 = %f\n',param.threshold2);
param.contour_int1 = centers(i);
param.contour_int2 = centers(j);
fprintf('contour_int1 = %f\n',param.contour_int1);
fprintf('contour_int2 = %f\n\n\n',param.contour_int2);

if param.plot
    f = figure;
    plot(centers,cdf_LFM1,'b');
    hold on;
    plot(centers,cdf_LFM2,'r');
    plot(centers(i),cdf_LFM1(i),'bo');
    plot(centers(j),cdf_LFM2(j),'ro');
    xlabel('Intensity of voxel [uint16]');
    ylabel('fraction of population');
    title('cumulative density function');
    legend('LFM1','LFM2');
    hold off;
    %arbitrary_scale = 4;
    %xlim([0 centers(arbitrary_scale*max(i,j))]);
    xlim([0 1e4]);
    mystr1 = [sprintf('%3.3f%% of voxels in LFM1 exceed %.0f in intensity\n',...
        100*(1-cdf_LFM1(i)), param.threshold1 )...
              sprintf('%3.3f%% of voxels in LFM2 exceed %.0f in intensity\n',...
        100*(1-cdf_LFM2(j)), param.threshold2 )...
        sprintf('\nLFM1 has %d voxels\n', numel(LFM1))...
        sprintf('LFM2 has %d voxels', numel(LFM2))...
        ];
   
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1);
    mystr2 = [num2str(param.N) ' bins in distribution'];
    text(0.4,0.3,mystr2,'FontSize',9,'Color',[0 0 0],'Interpreter','none');
    text(0.4,0.5,mystr1,'FontSize',9,'Color',[0 0 0],'Interpreter','none');
    
    str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
    save_plot(f,str);
end
end


function i = derive_threshold (cdf, param)
% look for large jump in cdf
dif = cdf(2:end)-cdf(1:end-1);
sdif = flipud(sortrows([dif' [1:length(dif)]'],1));
if sdif(1,1) > param.dynamic_range_thresh * sdif(2,1)
    i = sdif(1,2)+1;
else
    % if not found
    i = find(cdf>param.pop_thresh,1); % keep only first instance
end
end


function out = combineVols_iter_dim3 (LFM1, LFM2, param)
out = zeros(size(LFM2),'uint16');
% chunk LFM2 in third (z) dimension
% full frame in first (y) and second (x) dimensions
%M = 20; % number of chunks
s = size(LFM2);
%N = ceil(s(3)/M); % number of voxels/frames per chunk
%F = s(1)*s(2); % number of voxels per frame
fprintf('\n');
% for each slice of dim 3 in LFM2
parfor i = 1:s(3)
    %disp(sprintf('chunk %d of %d',i,M));
    % calculate index range
    zind = i;
    %zind = [(i-1)*N+1  i*N];
    if zind(2)>s(3)
        zind(2) = s(3);
    end
    linind = [(i-1)*N*F+1  i*N*F];
    if linind(1)>numel(LFM2)
        disp('LFM2 is complete!');
        continue;
    end
    if linind(2)>numel(LFM2)
        linind(2) = numel(LFM2);
    end
    fprintf('chunk %2d of %2d, z index %3d to %3d, linear index %12d to %12d\n',i,M, zind(1),zind(2),linind(1),linind(2));
    % for each voxel in chunk in LFM2, extract position to final_pos2
    %disp('Init_pos');
    final_pos2 = init_pos([linind(1):linind(2)]', LFM2, param);
    % for each position in final_pos2, apply rot and translation
    %disp('Translate');
    tmp = translate (final_pos2, -param.centroid);
    %disp('Rotate');
    tmp = rotate (tmp, param.rot); % CHECK
    %disp('Translate');
    final_pos2 = translate (tmp, param.trans); 
    % for each voxel in LFM1, find positions in final_pos2 and combine and normalize
    %disp('combine');
    out(i) = combine (LFM1, LFM2, final_pos2, linind, param);
end
if ~param.rapid
    outFile = sprintf('%s_%s.tif',param.timestamp,param.myfunc_combine);
%     if length(param.myfunc)>4 && strcmp(param.myfunc(end-3:end),'norm')
%         m = max(max(max(out)));
%         out = uint16( single(out)*2^16/single(m) );
%     end
    save_vol( out, param.savePath, outFile);
end

XguessSAVE1 = out;
outFile = sprintf('%s%s_%s.mat',param.savePath,param.timestamp,param.myfunc_combine);
fprintf('Saving combined volume to %s.\n',outFile);
save(outFile,'XguessSAVE1','-v7.3');
end



function out = combineVols (LFM1, LFM2, param)
out = zeros(size(LFM2),'uint16');
% chunk LFM2 in third (z) dimension
% full frame in first (y) and second (x) dimensions
M = 20; % number of chunks
s = size(LFM2);
N = ceil(s(3)/M); % number of voxels/frames per chunk
F = s(1)*s(2); % number of voxels per frame
fprintf('\n');
% for each chunk in LFM2
for i = 1:M
    %disp(sprintf('chunk %d of %d',i,M));
    % calculate index range
    zind = [(i-1)*N+1  i*N];
    if zind(2)>s(3)
        zind(2) = s(3);
    end
    linind = [(i-1)*N*F+1  i*N*F];
    if linind(1)>numel(LFM2)
        disp('LFM2 is complete!');
        break;
    end
    if linind(2)>numel(LFM2)
        linind(2) = numel(LFM2);
    end
    fprintf('chunk %2d of %2d, z index %3d to %3d, linear index %12d to %12d\n',i,M, zind(1),zind(2),linind(1),linind(2));
    % for each voxel in chunk in LFM2, extract position to final_pos2
    %disp('Init_pos');
    final_pos2 = init_pos([linind(1):linind(2)]', LFM2, param);
    % for each position in final_pos2, apply rot and translation
    %disp('Translate');
    tmp = translate (final_pos2, -param.centroid);
    %disp('Rotate');
    tmp = rotate (tmp, param.rot); % CHECK
    %disp('Translate');
    final_pos2 = translate (tmp, param.trans); 
    % for each voxel in LFM1, find positions in final_pos2 and combine and normalize
    %disp('combine');
    out(linind(1):linind(2)) = combine (LFM1, LFM2, final_pos2, linind, param);
end
if ~param.rapid
    outFile = sprintf('%s_%s.tif',param.timestamp,param.myfunc_combine);
%     if length(param.myfunc)>4 && strcmp(param.myfunc(end-3:end),'norm')
%         m = max(max(max(out)));
%         out = uint16( single(out)*2^16/single(m) );
%     end
    save_vol( out, param.savePath, outFile);
end

XguessSAVE1 = out;
outFile = sprintf('%s%s_%s.mat',param.savePath,param.timestamp,param.myfunc_combine);
fprintf('Saving combined volume to %s.\n',outFile);
save(outFile,'XguessSAVE1','-v7.3');
end

function out = loadData (f, param)
if exist(f,'file') == 2
    load(f);
    if exist('Xvolume','var')
        if isa(Xvolume,'uint16')
            XguessSAVE1 = Xvolume;
        else
            disp('Warning input data is not 16 bit.');
            keyboard
        end
        clear Xvolume;
    elseif exist('XguessSAVE1','var')
        disp('XguessSAVE1 found.');
        if ~isa(XguessSAVE1,'uint16')
            disp('Warning input data is not 16 bit.');
            if isa(XguessSAVE1,'uint8')
                disp('Data will be scaled to 16 bit.');
                XguessSAVE1 = cast(XguessSAVE1, 'uint16')*2^8;
            else
                keyboard
            end
        end
    else
        %disp('WTF?! Unknown data name.');
        fprintf('No data was recognized:\n\n');
        whos
        keyboard;
    end
    out = interpolate (XguessSAVE1,param);
else
    fprintf('File not found:\n%s\n\n',f);
    keyboard;
end 
end

function out = interpolate (LFM, param)
% initialize container of new size
s = size(LFM);
out = zeros(s(1),s(2),param.interp*s(3),'uint16');
boundary = param.interp/2 + 1;
for i=1:param.interp*s(3)
    if i < boundary
        out(:,:,i) = LFM(:,:,1);
    elseif i > (param.interp*s(3)-(boundary-1))
        out(:,:,i) = LFM(:,:,end);
    else
        a = round((i-1)/param.interp);
        b = a+1;
        N = 2*param.interp;
        fb = 1+2*mod(i-(param.interp-1),param.interp);
        fa = N-fb;
        out(:,:,i) = fa/N*LFM(:,:,a) + fb/N*LFM(:,:,b);
    end
end
end


function out = combine (LFM1, LFM2, chunk, linind, param)
% chunk = x,y,z centroids of voxels in LFM2
% for each centroid in chunk, calculate index of corresponding voxels in LFM1
L = length(chunk);
scale = [1/param.voxel_y*ones(L,1) ...
    1/param.voxel_x*ones(L,1) ...
    1/param.voxel_z*ones(L,1)];
abc_LFM1 = ceil(chunk.*scale);
% purge any indices that lie outLFM of LFM1
ind = find(abc_LFM1<1);
[r1,~] = ind2sub(size(abc_LFM1),ind);
s = size(LFM1);
ind = find(abc_LFM1(:,1)>s(1));
[r2,~] = ind2sub(size(abc_LFM1),ind);
ind = find(abc_LFM1(:,2)>s(2));
[r3,~] = ind2sub(size(abc_LFM1),ind);
ind = find(abc_LFM1(:,3)>s(3));
[r4,~] = ind2sub(size(abc_LFM1),ind);
rem = unique([r1;r2;r3;r4]);
N = length(abc_LFM1);
f_not_used = length(rem)/N;
fprintf('[Voxels] Number %d, fraction used %1.3f, with negative index %1.3f, outLFM of LFM1 %1.3f\n',...
    N,1-f_not_used,length(r1)/N,(length(r2)+length(r3)+length(r4))/N);
% chunk is a systematic sweep in z dimensino of LFM2
% in contrast, abc_LFM1 are the voxels in LFM1 coordinate system of a
% rotated and translated LFM2, so not well behaved
% thus, compute intensity in LFM2 coordinate system
out = ones(length(chunk),1,'uint16');
%out = zeros(length(chunk),1,'uint16');
%
% how to handle rem entries
% in abc, replace with 1,1,1
abc_LFM1(rem,:) = [ones(length(rem),1) ones(length(rem),1) ones(length(rem),1)];
%
% in intensity calculation, set all rem entries to zero
ind = sub2ind(size(LFM1),abc_LFM1(:,1),abc_LFM1(:,2),abc_LFM1(:,3));
i1 = LFM1(ind);
i1(rem) = 0.0;
i2 = LFM2(linind(1):linind(2))';

if strcmp(param.myfunc_combine,'multiply')
    % multiply
    out = uint16( single(i1).*single(i2) );
elseif strcmp(param.myfunc_combine,'multiply_sqrt')
    % multiply
    out = uint16( sqrt(single(i1).*single(i2)) );
elseif strcmp(param.myfunc_combine,'i1')
    out = i1;
elseif strcmp(param.myfunc_combine,'i2')
    out = i2;
elseif strcmp(param.myfunc_combine,'sum')
    out = i1 + i2;
elseif strcmp(param.myfunc_combine,'min')
    % OR min OR probability
    out = uint16( min(i1,i2));
else
    disp('WTF?!');
    keyboard
end
%max(max(max(out)))
end


function [cdf,centers,nullMIvec] = null_distribution (LFM1, LFM2, canonical, param)
nullf = [param.savePath param.inputFileName{1}(1:end-4) '_null.mat'];
if exist(nullf,'file') == 2
    load(nullf,'nullMIvec');
else
    nullMIvec = [];
end
LL = length(nullMIvec);
fprintf('\nnull distribution has N = %d\n',LL);
% determine gain so offset limit = half or quarter of volume limits
% and so that rotation limit = pi
a = max ([ size(LFM1) size(LFM2)]);
offset_limit = 0.15 * a * param.voxel_y;
gain = offset_limit / param.trans_amp;
tmp = param.rot_amp;
param.rot_amp = [pi/gain 0 0];
%N = param.Nnull;
%profile on;
fprintf('\nCount    Mutual_Information                Offset [um]                        Rotation [radians]\n');
parfor i=1:param.Nnull
    MI = 0;
    % perturb pos
    % randomly pick a translation vector and rotation vector
    % to be added to current location
    [d,r] = perturb(param,gain);
    % apply transformation
    % rotate an amount r PLUS param.rot
    % thus param.rot tracks the current rotation
    %rotated = rotate (canonical,r);
    % translate rotated by an amount param.trans+d
    % thus param.trans tracks the current position
    new = translate (rotate (canonical, r), param.centroid+d);
    %new = translate (rotated, param.centroid+d);
    % measure mutual_information
    if strcmp(param.myfunc_MI,'multiply')
        MI = mutual_information (LFM1, new, LFM2, param);
%     elseif strcmp(param.myfunc_MI,'multiply_sqrt')
%         MI = mutual_information_sqrt (LFM1, new, LFM2, param);
    else
        disp('WTF!');
        keyboard;
    end
    nullMIvec = [nullMIvec MI];
    str = sprintf('i = %4d, MI = %16.0f, d = [%7.3f0 %7.3f0 %7.3f0], r = [%7.3f0 %7.3f0 %7.3f0]',i,MI,d(1),d(2),d(3),r(1),r(2),r(3));
    disp(str);
    %         profile off
    %         profile viewer
end
param.rot_amp = tmp;

% add to existing null distribution if any
fprintf('\nnull distribution has N = %d, was N = %d\n',length(nullMIvec),LL);
save(nullf,'nullMIvec');

% CDF
centers = linspace(min(nullMIvec),max(nullMIvec),param.Nnull);
del = 0.5 * (centers(2)-centers(1));
edges  = [centers-del centers(end)+del];
%centers = 0.5 * (edges(1:end-1)+edges(2:end));
% this works. to check: sum(h_LFM1) == numel(LFM1)
h_MI = histcounts(nullMIvec,edges);
% normalize cdf to 1
h_MI = h_MI / sum(h_MI);

total = 0;
cdf = [];
for i=1:length(h_MI)
    total = total + h_MI(i);
    cdf = [cdf total];
end

% plot PDF and CDF
f = figure;
histogram(nullMIvec,50);
if strcmp(param.myfunc_MI,'multiply')
    xlabel('mutual information = sum(LFM1*LFM2)');
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    xlabel('mutual information = sum(sqrt(LFM1*LFM2))');
else
    disp('WTF!');
    keyboard;
end
ylabel('count');
title(['null distribution (bootstrapped from'...
    sprintf(' %d random registrations)',length(nullMIvec))]);

str=sprintf('%s%s_null_distribution.png',param.savePath,param.timestamp);
save_plot(f, str);

f = figure;
plot(centers,cdf);
if strcmp(param.myfunc_MI,'multiply')
    xlabel('mutual information = sum(LFM1*LFM2)');
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    xlabel('mutual information = sum(sqrt(LFM1*LFM2))');
else
    disp('WTF!');
    keyboard;
end
ylabel('normalized cumulative density function');
ylim([0 1]);
title(['null distribution (bootstrapped from'...
    sprintf(' %d random registrations)',length(nullMIvec))]);

str=sprintf('%s%s_null_cdf.png',param.savePath,param.timestamp);
save_plot(f, str);
end




function [new, param] = simulated_annealing (LFM1, new, LFM2, canonical, param)
%print_param(param);
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, new, LFM2, param);
% elseif strcmp(param.myfunc_MI,'multiply_sqrt')
%     MI = mutual_information_sqrt (LFM1, new, LFM2, param);
else
    disp('WTF!');
    keyboard;
end
last_MI = MI;
param = setT0 (MI,param);
param.MIvec = MI;
% set initial T. Start with T sufficiently high to "melt" the system
T = param.T0;
p = param.init_p;
param.Pvec = p;
% set max number of temperature changes and mean changes
Tchanges = param.TC0;
% while system not frozen and more temperature changes are allowed
% profile on;
param.transvec = [param.trans];
param.offsetvec = [param.offset];
param.rotvec = [param.rot];
% CDF
param.cdfvec = [];
w = find(param.centers>last_MI,1); % keep only first instance
if ~isempty(w)
    param.cdfvec = [param.cdfvec param.cdf(w)];
else
    val = double(last_MI)/double(max(param.nullMIvec));
    param.cdfvec = [param.cdfvec val] ;
end
disp('Count      Perturbation     Transformation      Overlap     Probability,Decision   ');
while Tchanges > 0
    pos_changes = param.MC0;
    while pos_changes > 0
        % perturb pos
        % randomly pick a translation vector and rotation vector
        % to be added to current location
        [d,r] = perturb(param,p);
        str0 = sprintf('d = [%7.3g0 %7.3g0 %7.3g0], r = [%7.3g0 %7.3g0 %7.3g0]',d(1),d(2),d(3),r(1),r(2),r(3));
        % apply transformation
        % rotate an amount r PLUS param.rot
        % thus param.rot tracks the current rotation
        rotated = rotate (canonical,param.rot+r);
        param.rot = param.rot + r;
        % translate rotated by an amount param.trans+d
        % thus param.trans tracks the current position
        new = translate (rotated, param.trans+d);
        param.trans = param.trans + d;
        str1 = print_param(param);
        % measure mutual_information
        if strcmp(param.myfunc_MI,'multiply')
            MI = mutual_information (LFM1, new, LFM2, param);
%         elseif strcmp(param.myfunc_MI,'multiply_sqrt')
%             MI = mutual_information_sqrt (LFM1, new, LFM2, param);
        else
            disp('WTF!');
            keyboard;
        end            
        delmi = MI - param.MIvec(end);
        str2 = sprintf('test MI = %7.3g, delmi = %7.3g',MI,delmi);
        if delmi > 0.0
            % keep
            str3 = 'Accept - MI increased';
            param.MIvec = [param.MIvec MI];
            %param.MItvec = [param.MItvec MIt];
            last_MI = MI;
        else
            %         # else accept D with probability P = exp(-E/kBT)
            %         # using (psuedo-)random number uniformly distributed in the interval (0,1)
            %p = calc_boltzman_p(MI,T,param);
            % Specifically, if random number is less the P, then accept
            rnd = rand(1);
            if rnd < p
                % accept decrease in MI mutual information
                str3 = sprintf('Accept %.3g < %.3g',rnd,p);
                param.MIvec = [param.MIvec MI];
                %param.MItvec = [param.MItvec MIt];
                last_MI = MI;
            else
                % reject move
                str3 = sprintf('Reject %.3g > %.3g',rnd,p);
                param.rot = param.rot - r;
                param.trans = param.trans - d;
                tmp = param.MIvec(end);
                param.MIvec = [param.MIvec tmp];
                %tmpt = param.MItvec(end);
                %param.MItvec = [param.MItvec tmpt];
                
            end
        end
        pos_changes = pos_changes-1;
        str4 = sprintf('%d %d',Tchanges,pos_changes);
        str6 = sprintf('Final MI = %.3g',last_MI);
        str5 = print_param(param);
        val7 = param.trans - param.centroid;
        str7 = sprintf('final offset = [%6.6g0 %6.6g0 %6.6g0]',val7(1),val7(2),val7(3));
        w = find(param.centers>last_MI,1); % keep only first instance
        if last_MI > max(param.nullMIvec)
            val = double(last_MI)/double(max(param.nullMIvec));
            str8 = sprintf('MI frac = %.5g, siman exceeds null',val);
            param.cdfvec = [param.cdfvec val] ;
        elseif ~isempty(w)
            str8 = sprintf('MI frac = %.5g, null exceeds siman',param.cdf(w));
            param.cdfvec = [param.cdfvec param.cdf(w)];
        else
            disp('WTF?');
            keyboard
        end
        dif = last_MI - max(param.nullMIvec);
        % last_MI > max(param.nullMIvec) =>dif>0
        % last_MI < max(param.nullMIvec) =>dif<0
%         if dif > 0
%             str9 = sprintf('siman exceeds null by %3.2f%%',100*abs(dif)/max(param.nullMIvec));
%         else
%             str9 = sprintf('null exceeds siman by %3.2f%%',100*abs(dif)/last_MI);
%         end
        fprintf('%7s%75s  %84s  %40s  %22s  %84s  %20s %40s %20s\n',str4,str0,str1,str2,str3,str5,str6,str7,str8);
        %fprintf('%7s%75s  %84s  %40s  %22s  %84s  %20s %40s %20s %20s\n',str4,str0,str1,str2,str3,str5,str6,str7,str8,str9);
        param.Pvec = [param.Pvec p];
        param.transvec = [param.transvec; param.trans];
        param.offsetvec = [param.offsetvec; val7];
        param.rotvec = [param.rotvec; param.rot];
%         profile off
%          profile viewer
%          keyboard
    end
    %profile viewer;
    %keyboard
    T = lowerT(T,param);
    Tchanges = Tchanges-1;
    p = p * param.prate;
end
rotated = rotate (canonical,param.rot);
new = translate (rotated, param.trans);
end

function save_plot (f, filename)
print(f,filename,'-dpng');
end

function save_stats (param)
prefix = sprintf('%s%s_',param.savePath,param.timestamp);
% plot Pvec
h = figure;
plot(1:numel(param.Pvec),param.Pvec);
xlabel('iteration');
ylabel('probability of allowing a decrease in MI');
title('simulated annealing schedule');
if 0>1
    fname = sprintf('%s_p.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_p.png',prefix);
    save_plot(h, str);
end

h = figure;
plot(1:numel(param.MIvec),log10(param.MIvec));
xlabel('iteration');
ylabel('log (mutual information)');
title('evolution of mutual information during simulated annealing');
if 0>1
    fname = sprintf('%s_MI.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_MI.png',prefix);
    save_plot(h, str);
end

f = figure;
plot(1:numel(param.cdfvec),param.cdfvec);
hold on;
plot(xlim,[1 1],'--r');
hold off;
xlabel('iteration');
ylabel('fraction of null distribution');
title('fraction of null distribution less than current mutual informaton');
if 0>1
    fname = sprintf('%s_frac_null.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_frac_null.png',prefix);
    save_plot(f, str);
end

h = figure;
a = size(param.offsetvec);
subplot(1,3,1);
plot(1:a(1),param.offsetvec(:,1));
xlabel('iteration');
ylabel('offset in dim 1 (um)');
hold on;
plot(1,param.offsetvec(1,1),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,2);
plot(1:a(1),param.offsetvec(:,2));
xlabel('iteration');
ylabel('offset in dim 2 (um)');
hold on;
plot(1,param.offsetvec(1,2),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,3);
plot(1:a(1),param.offsetvec(:,3));
xlabel('iteration');
ylabel('offset in dim 3 (um)');
hold on;
plot(1,param.offsetvec(1,3),'ro');
hold off;
xlim([1 a(1)]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'trajectory of simulated annealing';
text(0.4,0.97,str,'FontSize',12,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_offset.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_offset.png',prefix);
    save_plot(h, str);
end

h = figure;
a = size(param.rotvec);
subplot(1,3,1);
plot(1:a(1),param.rotvec(:,1));
xlabel('iteration');
ylabel('rotation around dim 1 (radians)');
hold on;
plot(1,param.rotvec(1,1),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,2);
plot(1:a(1),param.rotvec(:,2));
xlabel('iteration');
ylabel('rotation around dim 2 (radians)');
%title('trajectory of simulated annealing');
hold on;
plot(1,param.rotvec(1,2),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,3);
plot(1:a(1),param.rotvec(:,3));
xlabel('iteration');
ylabel('rotation around dim 3 (radians)');
hold on;
plot(1,param.rotvec(1,3),'ro');
hold off;
xlim([1 a(1)]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'trajectory of simulated annealing';
text(0.4,0.98,str,'FontSize',12,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_rot.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_rot.png',prefix);
    save_plot(h, str);
end

end


function out = print_param (param)
a=param.trans;
b=param.rot;
out = sprintf('trans = [%6.6g0 %6.6g0 %6.6g0], rot = [%7.6g0 %7.6g0 %7.6g0]',a(1),a(2),a(3),b(1),b(2),b(3));
end


function mi = mutual_information (LFM1, new, LFM2, param)
scale = [1/param.voxel_y 0 0 ; 0 1/param.voxel_x 0 ; 0 0 1/param.voxel_z];
new_pixels = ceil(new*scale);
s = size(LFM1);
% find rows of new_pixels that overlap LFM1 and save as index
index = find( new_pixels(:,1)<=s(1) & new_pixels(:,2)<=s(2) & new_pixels(:,3)<=s(3) ...
    & new_pixels(:,1)>0 & new_pixels(:,2)>0 & new_pixels(:,3)>0 );
% in LFM1, lookup intensity value at 
i1 = LFM1(sub2ind(s,new_pixels(index,1),new_pixels(index,2),new_pixels(index,3)));
i2 = LFM2(param.index2(index));
mi = sum(double(i1).*double(i2));
end


function out = rotation_matrix (angle)
a = [1 0 0;...
    0 cos(angle(1)) -sin(angle(1));...
    0 sin(angle(1)) cos(angle(1))];
b = [cos(angle(2)) 0 -sin(angle(2));...
    0 1 0;...
    sin(angle(2)) 0 cos(angle(2))];
c = [cos(angle(3)) -sin(angle(3)) 0;...
    sin(angle(3)) cos(angle(3)) 0;...
    0 0 1];
out = a*b*c;
end


function [out] = rotate (pos, angle)
rot = single( rotation_matrix (angle) );
out = pos*rot;
end


function [out] = translate (pos, delta)
D = ones(size(pos),'single').*delta;
out = pos + D;
end

% function p = calc_boltzman_p (T, param)
% val = -param.scale/T;
% if val > 0.0
%     p = 1.0;
% else
%     p = exp(val);
% end
% if p<0.0 || p>1.0
%     disp('WTF?');
%     keyboard
% end
% end


function out = calc_centroid (LFM, param)
s = size(LFM);
out = s/2.*[param.voxel_y param.voxel_x param.voxel_z];
end


function pos = init_pos (linear_i,LFM, param)
[a,b,c] = ind2sub(size(LFM),linear_i);
y = single( (a-0.5) * param.voxel_y );
x = single( (b-0.5) * param.voxel_x );
z = single( (c-0.5) * param.voxel_z );
pos = [y x z];
end



function save_vol (A, savePath, outFile)
imwrite( squeeze(A(:,:,1)), [savePath outFile]);
for k = 2:size(A,3)
    imwrite(squeeze(A(:,:,k)),  [savePath outFile], 'WriteMode', 'append');
end
end


function [d,r] = perturb (param, gain)
d1 = random('unif',-1,1) * param.trans_amp * gain;
d2 = random('unif',-1,1) * param.trans_amp * gain;
d3 = random('unif',-1,1) * param.trans_amp * gain;
r1 = random('unif',-1,1) * param.rot_amp(1) * gain;
r2 = random('unif',-1,1) * param.rot_amp(2) * gain;
r3 = random('unif',-1,1) * param.rot_amp(3) * gain;
d = [d1 d2 d3];
r = [r1 r2 r3];
%r = [r1 0 0];
end


function out =  lowerT(T,param)
delT = -param.Trate*T;
out = T+delT;
end


function param = setT0 (MI, param)
% aim for 0.95 aceptance at T0
param.T0 = -1/log(param.init_p);
%param.TC0 = round(log10(param.final_p) / log10(param.Trate))
%param.T0 = -1/log(0.95)/MI;
param.prate = 10^( ( log10(param.final_p)-log10(param.init_p) ) / param.TC0);
end



function print_fraction (index, LFM, str)
i = 0;
for j=1:length(index)
    i = i + single(LFM(index(j)));
end
iT = single(sum(sum(sum(LFM))));

fi = i/iT;
fv = single(length(index))/single(numel(LFM));
fprintf('%s: %d of %d voxels (%7.3g), %d of %d total intensity (%7.3g)\n',...
    str,length(index),numel(LFM),fv,...
    i,iT,fi);
end


function offsets = estimate_offsets (LFM1, LFM2, new, param)

onetwo = squeeze(max(LFM1,[],3));
twothree = squeeze(max(LFM1,[],1));
LFM1_d1 = single(squeeze(max(onetwo,[],2)));
LFM1_d2 = single(squeeze(max(onetwo,[],1)));
LFM1_d3 = single(squeeze(max(twothree,[],1)));

%
% projections for new
%

% convert new from microns to pixels
a = size(new);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_new = ceil(new.*scale);
ns = max(abc_new);

% create 2d projections of new in pixel space
onetwo = zeros(ns(1),ns(2));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,2) > 0
        onetwo(abc_new(i,1),abc_new(i,2)) = max([onetwo(abc_new(i,1),abc_new(i,2)) double(LFM2(param.index2(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(new)
    if abc_new(i,2) > 0 && abc_new(i,3) > 0
        twothree(abc_new(i,2),abc_new(i,3)) = max([twothree(abc_new(i,2),abc_new(i,3))  double(LFM2(param.index2(i)))]);
    end
end
% convert 2d projections to 1d
new_d1 = squeeze(max(onetwo,[],2));
new_d2 = squeeze(max(onetwo,[],1));
new_d3 = squeeze(max(twothree,[],1));

%
% xcorr
%

% [u_acor, u_lag]  = xcorr( LFM1_d1,   new_d1 );
% [v_acor, v_lag]  = xcorr( LFM1_d2,   new_d2 );
% [w_acor, w_lag]  = xcorr( LFM1_d3,   new_d3 );
% 
% u_index = find(u_acor == max(u_acor));
% v_index = find(v_acor == max(v_acor));
% w_index = find(w_acor == max(w_acor));
% 
% offsets = [u_lag(u_index)*param.voxel_y...
%     v_lag(v_index)*param.voxel_x...
%     w_lag(w_index)*param.voxel_z...
%     ];
% 
% if param.plot
%     f = figure;
%     subplot(3,1,1);
%     plot(u_lag,u_acor/max(u_acor));
%     xlabel('lag in first dimension [pixels]');
%     ylabel('correlation');
%     title('Normalized cross-correlation between LFM1 and rotated LFM2');
%     subplot(3,1,2);
%     plot(v_lag,v_acor/max(v_acor));
%     xlabel('lag in second dimension [pixels]');
%     ylabel('correlation');
%     subplot(3,1,3);
%     plot(w_lag,w_acor/max(w_acor));
%     xlabel('lag in third dimension [pixels]');
%     ylabel('correlation');
%     
%     str=sprintf('%s%s_xcorr.png',param.savePath,param.timestamp);
%     save_plot(f, str);
% end

LFM1_d1_log10 = log10(LFM1_d1);
LFM1_d1_log10(find(isinf(LFM1_d1_log10))) = 0;
LFM1_d2_log10 = log10(LFM1_d2);
LFM1_d2_log10(find(isinf(LFM1_d2_log10))) = 0;
LFM1_d3_log10 = log10(LFM1_d3);
LFM1_d3_log10(find(isinf(LFM1_d3_log10))) = 0;

new_d1_log10 = log10(new_d1);
new_d1_log10(find(isinf(new_d1_log10))) = 0;
new_d2_log10 = log10(new_d2);
new_d2_log10(find(isinf(new_d2_log10))) = 0;
new_d3_log10 = log10(new_d3);
new_d3_log10(find(isinf(new_d3_log10))) = 0;

[u_acor, u_lag]  = xcorr( LFM1_d1_log10,   new_d1_log10 );
[v_acor, v_lag]  = xcorr( LFM1_d2_log10,   new_d2_log10 );
[w_acor, w_lag]  = xcorr( LFM1_d3_log10,   new_d3_log10 );

u_index = find(u_acor == max(u_acor));
v_index = find(v_acor == max(v_acor));
w_index = find(w_acor == max(w_acor));

offsets = [u_lag(u_index)*param.voxel_y...
    v_lag(v_index)*param.voxel_x...
    w_lag(w_index)*param.voxel_z...
    ];

if param.plot
    f = figure;
    subplot(3,1,1);
    plot(u_lag,u_acor/max(u_acor));
    hold on;
    plot([u_lag(u_index) u_lag(u_index)],[0 1],'r');
    hold off;
    xlabel('lag in first dimension [pixels]');
    ylabel('correlation');
    title('Normalized log cross-correlation between LFM1 and rotated LFM2');
    str = sprintf('offset = %.0d pixels %.1f um',u_lag(u_index),u_lag(u_index)*param.voxel_y);
    text(100,0.2,str,'FontSize',8,'Color',[0 0 0] ,'Interpreter','none');
    ylim([0 1]);
    subplot(3,1,2);
    plot(v_lag,v_acor/max(v_acor));
    hold on;
    plot([v_lag(v_index) v_lag(v_index)],ylim,'r');
    hold off;
    xlabel('lag in second dimension [pixels]');
    ylabel('correlation');
    str = sprintf('offset = %.0d pixels %.1f um',v_lag(v_index),v_lag(v_index)*param.voxel_x);
    text(100,0.2,str,'FontSize',8,'Color',[0 0 0] ,'Interpreter','none');
    ylim([0 1]);
    subplot(3,1,3);
    plot(w_lag,w_acor/max(w_acor));
    hold on;
    plot([w_lag(w_index) w_lag(w_index)],ylim,'r');
    hold off;
    xlabel('lag in third dimension [pixels]');
    ylabel('correlation');
    str = sprintf('offset = %.0d pixels %.1f um',w_lag(w_index),w_lag(w_index)*param.voxel_z);
    text(100,0.2,str,'FontSize',8,'Color',[0 0 0] ,'Interpreter','none');
    ylim([0 1]);
    
    str=sprintf('%s%s_xcorr.png',param.savePath,param.timestamp);
    save_plot(f, str);
end

end



function param = save_1d_max_projections (LFM1, LFM2, pos1, pos2, new, param, str)
colors = {[0 0 1]
    [0.8 0.8 0.8]
    [0 0 1]
    [1 0 0]};

onetwo = squeeze(max(LFM1,[],3));
twothree = squeeze(max(LFM1,[],1));
LFM1_d1 = single(squeeze(max(onetwo,[],2)));
LFM1_d2 = single(squeeze(max(onetwo,[],1)));
LFM1_d3 = single(squeeze(max(twothree,[],1)));

onetwo = squeeze(max(LFM2,[],3));
twothree = squeeze(max(LFM2,[],1));
LFM2_d1 = single(squeeze(max(onetwo,[],2)));
LFM2_d2 = single(squeeze(max(onetwo,[],1)));
LFM2_d3 = single(squeeze(max(twothree,[],1)));

% projections for new
% convert new from microns to pixels
a = size(new);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_new = ceil(new.*scale);
ns = max(abc_new);

% create 2d projections of new in pixel space
onetwo = zeros(ns(1),ns(2));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,2) > 0
        onetwo(abc_new(i,1),abc_new(i,2)) = max([onetwo(abc_new(i,1),abc_new(i,2)) double(LFM2(param.index2(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(new)
    if abc_new(i,2) > 0 && abc_new(i,3) > 0
        twothree(abc_new(i,2),abc_new(i,3)) = max([twothree(abc_new(i,2),abc_new(i,3))  double(LFM2(param.index2(i)))]);
    end
end
% convert 2d projections to 1d
new_d1 = squeeze(max(onetwo,[],2));
new_d2 = squeeze(max(onetwo,[],1));
new_d3 = squeeze(max(twothree,[],1));

% projections for pos1
a = size(pos1);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_pos1 = ceil(pos1.*scale);
ns = max(abc_pos1);

onetwo = zeros(ns(1),ns(2));
for i=1:length(pos1)
    if abc_pos1(i,1) > 0 && abc_pos1(i,2) > 0
        onetwo(abc_pos1(i,1),abc_pos1(i,2)) = max([onetwo(abc_pos1(i,1),abc_pos1(i,2)) double(LFM1(param.index1(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(pos1)
    if abc_pos1(i,2) > 0 && abc_pos1(i,3) > 0
        twothree(abc_pos1(i,2),abc_pos1(i,3)) = max([twothree(abc_pos1(i,2),abc_pos1(i,3)) double(LFM1(param.index1(i)))]);
    end
end
pos1_d1 = squeeze(max(onetwo,[],2));
pos1_d2 = squeeze(max(onetwo,[],1));
pos1_d3 = squeeze(max(twothree,[],1));

% projections for pos2
a = size(pos2);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_pos2 = ceil(pos2.*scale);
ns = max(abc_pos2);

onetwo = zeros(ns(1),ns(2));
for i=1:length(pos2)
    if abc_pos2(i,1) > 0 && abc_pos2(i,2) > 0
        onetwo(abc_pos2(i,1),abc_pos2(i,2)) = max([onetwo(abc_pos2(i,1),abc_pos2(i,2)) double(LFM2(param.index2(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(pos2)
    if abc_pos2(i,2) > 0 && abc_pos2(i,3) > 0
        twothree(abc_pos2(i,2),abc_pos2(i,3)) = max([twothree(abc_pos2(i,2),abc_pos2(i,3)) double(LFM2(param.index2(i)))]);
    end
end
pos2_d1 = squeeze(max(onetwo,[],2));
pos2_d2 = squeeze(max(onetwo,[],1));
pos2_d3 = squeeze(max(twothree,[],1));

if 0 < 1
    T = '      (scaled by total intensity in sample)';
    scale_d1 = [sum(pos1_d1) sum(pos2_d1) sum(LFM1_d1) sum(LFM2_d1) sum(new_d1)];
    scale_d2 = [sum(pos1_d2) sum(pos2_d2) sum(LFM1_d2) sum(LFM2_d2) sum(new_d2)];
    scale_d3 = [sum(pos1_d3) sum(pos2_d3) sum(LFM1_d3) sum(LFM2_d3) sum(new_d3)];
    yl = 'normalized intensity';
    %disp(T);
elseif 0 < 1
    T = '      (scaled by max intensity in each projected volume)';
    scale_d1 = single([max(pos1_d1) max(pos2_d1) max(LFM1_d1) max(LFM2_d1) max(new_d1)]);
    scale_d2 = single([max(pos1_d2) max(pos2_d2) max(LFM1_d2) max(LFM2_d2) max(new_d2)]);
    scale_d3 = single([max(pos1_d3) max(pos2_d3) max(LFM1_d3) max(LFM2_d3) max(new_d3)]);
    yl = 'normalized intensity';
    %disp(T);
else
    T = '      (not scaled)';
    scale_d1 = [1 1 1 1 1];
    scale_d2 = [1 1 1 1 1];
    scale_d3 = [1 1 1 1 1];
    yl = 'intensity';
    %disp(T);
end

f = figure;
set(gcf,'Position',[79          18        1270         940]);
subplot(3,1,1);
semilogy(LFM1_d1/scale_d1(3),'.','Color',colors{1});
hold on;
semilogy(LFM2_d1/scale_d1(4),'.','Color',colors{2});
semilogy(new_d1/scale_d1(5),'-','Color',colors{4});
xlabel('first dimension [pixels]');
ylabel(yl);
hold off;
legend('LFM1','LFM2','LFM2 coarse reg');
title([str T],'Interpreter','none');
if param.xlim_1d(1) > 0
    xlim([0 param.xlim_1d(1)]);
else
    a = xlim;
    xlim([0 a(2)]);
    param.xlim_1d(1)=a(2);
end

subplot(3,1,2);
semilogy(LFM1_d2/scale_d2(3),'.','Color',colors{1});
hold on;
semilogy(LFM2_d2/scale_d2(4),'.','Color',colors{2});
semilogy(new_d2/scale_d2(5),'-','Color',colors{4});
xlabel('second dimension [pixels]');
ylabel(yl);
hold off;
legend('LFM1','LFM2','LFM2 coarse reg');
if param.xlim_1d(2) > 0
    xlim([0 param.xlim_1d(2)]);
else
    a = xlim;
    xlim([0 a(2)]);
    param.xlim_1d(2)=a(2);
end

subplot(3,1,3);
semilogy(LFM1_d3/scale_d3(3),'.','Color',colors{1});
hold on;
semilogy(LFM2_d3/scale_d3(4),'.','Color',colors{2});
semilogy(new_d3/scale_d3(5),'-','Color',colors{4});
xlabel('third dimension [pixels]');
ylabel(yl);
hold off;
legend('LFM1','LFM2','LFM2 coarse reg');
if param.xlim_1d(3) > 0
    xlim([0 param.xlim_1d(3)]);
else
    a = xlim;
    xlim([0 a(2)]);
    param.xlim_1d(3)=a(2);
end

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
save_plot(f, str);
end


function save_2d_contour_plots (LFM1, LFM2, pos1, pos2, new, param, str)
colors = {[1 0 0]
    [0.8 0.8 0.8]
    [0 0 1]
    [0.4 0.4 0.4]};

% projections for new
% convert new from microns to pixels
a = size(new);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_new = ceil(new.*scale);
ns = max(abc_new);

% create 2d projections of new in pixel space
onetwo_new = zeros(ns(1),ns(2));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,2) > 0
        onetwo_new(abc_new(i,1),abc_new(i,2)) ...
            = max([onetwo_new(abc_new(i,1),abc_new(i,2)) ...
            double(LFM2(param.index2(i)))]);
    end
end
onethree_new = zeros(ns(1),ns(3));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,3) > 0
        onethree_new(abc_new(i,1),abc_new(i,3)) ...
            = max([onethree_new(abc_new(i,1),abc_new(i,3))...
            double(LFM2(param.index2(i)))]);
    end
end
threetwo_new = zeros(ns(3),ns(2));
for i=1:length(new)
    if abc_new(i,3) > 0 && abc_new(i,2) > 0
        threetwo_new(abc_new(i,3),abc_new(i,2))...
            = max([threetwo_new(abc_new(i,3),abc_new(i,2))...
            double(LFM2(param.index2(i)))]);
    end
end


% projections for pos1
% convert new from microns to pixels
a = size(pos1);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_pos1 = ceil(pos1.*scale);
ns = max(abc_pos1);

% create 2d projections of new in pixel space
onetwo_pos1 = zeros(ns(1),ns(2));
for i=1:length(pos1)
    if abc_pos1(i,1) > 0 && abc_pos1(i,2) > 0
        onetwo_pos1(abc_pos1(i,1),abc_pos1(i,2))...
            = max([onetwo_pos1(abc_pos1(i,1),abc_pos1(i,2))...
            double(LFM1(param.index1(i)))]);
    end
end
onethree_pos1 = zeros(ns(1),ns(3));
for i=1:length(pos1)
    if abc_pos1(i,1) > 0 && abc_pos1(i,3) > 0
        onethree_pos1(abc_pos1(i,1),abc_pos1(i,3))...
            = max([onethree_pos1(abc_pos1(i,1),abc_pos1(i,3))...
            double(LFM1(param.index1(i)))]);
    end
end
threetwo_pos1 = zeros(ns(3),ns(2));
for i=1:length(pos1)
    if abc_pos1(i,3) > 0 && abc_pos1(i,2) > 0
        threetwo_pos1(abc_pos1(i,3),abc_pos1(i,2))...
            = max([threetwo_pos1(abc_pos1(i,3),abc_pos1(i,2))...
            double(LFM1(param.index1(i)))]);
    end
end

% projections for pos2
% convert new from microns to pixels
a = size(pos2);
scale = [1/param.voxel_y*ones(a(1),1) ...
    1/param.voxel_x*ones(a(1),1) ...
    1/param.voxel_z*ones(a(1),1)];
abc_pos2 = ceil(pos2.*scale);
ns = max(abc_pos2);

% create 2d projections of new in pixel space
onetwo_pos2 = zeros(ns(1),ns(2));
for i=1:length(pos2)
    if abc_pos2(i,1) > 0 && abc_pos2(i,2) > 0
        onetwo_pos2(abc_pos2(i,1),abc_pos2(i,2)) ...
            = max([onetwo_pos2(abc_pos2(i,1),abc_pos2(i,2)) ...
            double(LFM2(param.index2(i)))]);
    end
end
onethree_pos2 = zeros(ns(1),ns(3));
for i=1:length(pos2)
    if abc_pos2(i,2) > 0 && abc_pos2(i,3) > 0
        onethree_pos2(abc_pos2(i,1),abc_pos2(i,3)) ...
            = max([onethree_pos2(abc_pos2(i,1),abc_pos2(i,3)) ...
            double(LFM2(param.index2(i)))]);
    end
end
threetwo_pos2 = zeros(ns(3),ns(2));
for i=1:length(pos2)
    if abc_pos2(i,3) > 0 && abc_pos2(i,2) > 0
        threetwo_pos2(abc_pos2(i,3),abc_pos2(i,2)) ...
            = max([threetwo_pos2(abc_pos2(i,3),abc_pos2(i,2)) ...
            double(LFM2(param.index2(i)))]);
    end
end

% centroid in pixels
c = [param.centroid(1)/param.voxel_y...
    param.centroid(2)/param.voxel_x...
    param.centroid(3)/param.voxel_z];

f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
hold on;
a = size(onethree_new);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onethree_new,v,'Color',colors{1});
a = size(onethree_pos1);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,onethree_pos1,v,'Color',colors{3});
a = size(onethree_pos2);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onethree_pos2,v,'Color',colors{2});
xlabel('dim three [pixels]');
ylabel('dim one [pixels]');
plot(c(3),c(1),'g*')
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');
set(gca,'XDir','reverse');


subplot(2,2,2);
hold on;
a = size(onetwo_new);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onetwo_new,v,'Color',colors{1});
a = size(onetwo_pos1);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,onetwo_pos1,v,'Color',colors{3});
a = size(onetwo_pos2);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onetwo_pos2,v,'Color',colors{2});
xlabel('dim two [pixels]');
ylabel('dim one [pixels]');
plot(c(2),c(1),'g*')
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

subplot(2,2,4);
hold on;
a = size(threetwo_new);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,threetwo_new,v,'Color',colors{1});
a = size(threetwo_pos1);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,threetwo_pos1,v,'Color',colors{3});
a = size(threetwo_pos2);
x = 1:a(2);
y = 1:a(1);
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,threetwo_pos2,v,'Color',colors{2});
xlabel('dim two [pixels]');
ylabel('dim three [pixels]');
plot(c(2),c(3),'g*');
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);

text(0.2,0.45,str,'FontSize',12,'Color',[0 0 0] ,'Interpreter','none');

msg = sprintf('LFM1: contour at intensity level %4.0f',param.contour_int1);
text(0.2,0.4,msg,'FontSize',12,'Color',colors{3},'Interpreter','none');

msg = sprintf('LFM2: contour at intensity level %4.0f',param.contour_int2);
text(0.2,0.35,msg,'FontSize',12,'Color',colors{2},'Interpreter','none');

msg = sprintf('LFM2 coarse reg: contour at intensity level %4.0f',param.contour_int2);
text(0.2,0.3,msg,'FontSize',12,'Color',colors{1},'Interpreter','none');

text(0.2,0.25,'centroid of LFM2 coarse reg','FontSize',12,'Color',[0 1 0],'Interpreter','none');

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
save_plot(f, str);
end


function save_2d_max_projections (LFM1, LFM2, new, param, flag, str)

f = figure('units','normalized','outerposition',[0 0 1 1]);
yz = squeeze(max(LFM1,[],2));
subplot(2,6,1);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM1');
%colorbar();
daspect([1,1,1]);
set(gca,'XDir','reverse');

xy = squeeze(max(LFM1,[],3));
subplot(2,6,2);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM1');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(LFM1,[],1))';
subplot(2,6,8);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
title('LFM1');
colorbar();
daspect([1,1,1]);

yz = squeeze(max(LFM2,[],2));
subplot(2,6,3);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM2');
%colorbar();
daspect([1,1,1]);
set(gca,'XDir','reverse');


xy = squeeze(max(LFM2,[],3));
subplot(2,6,4);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM2');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(LFM2,[],1))';
subplot(2,6,10);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
title('LFM2');
colorbar();
daspect([1,1,1]);

if flag
    yz = squeeze(max(new,[],2));
    subplot(2,6,5);
    dr = ceil(log2(single(max(max(yz)))));
    imagesc(yz,[0 2^dr]);
    xlabel('three [pixels]');
    ylabel('one [pixels]');
    title('DLFM');
    %colorbar();
    daspect([1,1,1]);
    set(gca,'XDir','reverse');
    
    xy = squeeze(max(new,[],3));
    subplot(2,6,6);
    dr = ceil(log2(single(max(max(xy)))));
    imagesc(xy,[0 2^dr]);
    xlabel('two [pixels]');
    ylabel('one [pixels]');
    title('DLFM');
    %colorbar();
    daspect([1,1,1]);
    
    xz = squeeze(max(new,[],1))';
    subplot(2,6,12);
    dr = ceil(log2(single(max(max(xz)))));
    imagesc(xz,[0 2^dr]);
    xlabel('two [pixels]');
    ylabel('three [pixels]');
    title('DLFM');
    colorbar();
    daspect([1,1,1]);
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.1,0.98,['LFM1 = ' param.inputFilePath1 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
text(0.1,0.96,['LFM2 = ' param.inputFilePath2 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
if flag
    text(0.1,0.94,['DLFM = ' param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
end
text(0.4,0.1,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
save_plot(f, str);

end


function save_2d_max_projections_compact (LFM1, LFM2, new, param, flag, str)

f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1);
yz = squeeze(max(LFM1,[],2));
xy = squeeze(max(LFM1,[],3));
xz = squeeze(max(LFM1,[],1))';
a = size(LFM1);
big_image = [a(3)+a(1) a(3)+a(2)];
imagesc(zeros(big_image));
a = size(yz);
x = [1 a(2)];
y = [1 a(1)];
imagesc('XData',x,'YData',y,'CData',fliplr(yz));
b = size(xy);
x = [a(2)+1 b(2)+a(2)];
y = [1 b(1)];
imagesc('XData',x,'YData',y,'CData',xy);
c = size(xz);
x = [a(2)+1 b(2)+a(2)];
y = [b(1)+1 b(1)+c(1)];
imagesc('XData',x,'YData',y,'CData',xz);
title('LFM1');
colorbar();
xlabel('pixels');
ylabel('pixels');
daspect([1,1,1]);

subplot(1,3,2);
yz = squeeze(max(LFM2,[],2));
xy = squeeze(max(LFM2,[],3));
xz = squeeze(max(LFM2,[],1))';
a = size(LFM2);
big_image = [a(3)+a(1) a(3)+a(2)];
imagesc(zeros(big_image));
a = size(yz);
x = [1 a(2)];
y = [1 a(1)];
imagesc('XData',x,'YData',y,'CData',fliplr(yz));
b = size(xy);
x = [a(2)+1 b(2)+a(2)];
y = [1 b(1)];
imagesc('XData',x,'YData',y,'CData',xy);
c = size(xz);
x = [a(2)+1 b(2)+a(2)];
y = [b(1)+1 b(1)+c(1)];
imagesc('XData',x,'YData',y,'CData',xz);
title('LFM2');
colorbar();
xlabel('pixels');
ylabel('pixels');
daspect([1,1,1]);

if flag
    subplot(1,3,3);
    yz = squeeze(max(new,[],2));
    xy = squeeze(max(new,[],3));
    xz = squeeze(max(new,[],1))';
    a = size(new);
    big_image = [a(3)+a(1) a(3)+a(2)];
    imagesc(zeros(big_image));
    a = size(yz);
    x = [1 a(2)];
    y = [1 a(1)];
    imagesc('XData',x,'YData',y,'CData',fliplr(yz));
    b = size(xy);
    x = [a(2)+1 b(2)+a(2)];
    y = [1 b(1)];
    imagesc('XData',x,'YData',y,'CData',xy);
    c = size(xz);
    x = [a(2)+1 b(2)+a(2)];
    y = [b(1)+1 b(1)+c(1)];
    imagesc('XData',x,'YData',y,'CData',xz);
    title('DLFM');
    colorbar();
    xlabel('pixels');
    ylabel('pixels');
    daspect([1,1,1]);
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.1,0.98,['LFM1 = ' param.inputFilePath1 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
text(0.1,0.96,['LFM2 = ' param.inputFilePath2 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
if flag
    text(0.1,0.94,['DLFM = ' param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
end
text(0.4,0.1,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
save_plot(f, str);

end
