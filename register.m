function MI = register (param)
MI = -1;

%%
if ~param.justCalcMI
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
end

param

%% load volumes
f = [param.inputFilePath1 param.inputFileName{1}];
fprintf('\nLoading volume %s\n\n',f);
LFM1 = loadData(f, param);
f = [param.inputFilePath2 param.inputFileName{1}];
fprintf('\nLoading volume %s\n\n',f);
LFM2 = loadData(f, param);
param.voxel_z = param.voxel_z / param.interp;

%keyboard

%% if just transforming confocal data
if param.rapid && param.confocal
    tmp = LFM1;
    LFM1 = LFM2;
    LFM2 = ones(size(LFM2),'uint16');
    param.rot = -param.rot;
    param.centroid = calc_centroid(LFM2,param);
    param.trans = param.centroid;
end

%% if registration parameters are already known
%  then combine registered volumes
if param.rapid
    % must specify: param.centroid, param.trans, param.rot
    combineVols_iter_dim3 (LFM1, LFM2, param);
    %MI = mutual_information (LFM1, out, LFM2, param)
    param
    elapsedTime = toc
    diary off;
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

%keyboard
% fprintf('Total 

%% calculate thresholds
param = calculate_thresholds (LFM1, LFM2, param,'cdf_voxel_intensities');

%keyboard

%% voxel positions
param.index1 = find(LFM1>param.threshold1);
print_fraction(param.index1,LFM1,'LFM1');
pos1 = init_pos(param.index1,LFM1,param);
param.index2 = find(LFM2>param.threshold2);
print_fraction(param.index2,LFM2,'LFM2');
%param.size = size(LFM2);
pos2 = init_pos(param.index2,LFM2,param);

%% if justCalcMI
if param.justCalcMI
    param.centroid = calc_centroid(LFM2,param);
    tmp = translate (pos2, -param.centroid);
    tmp = rotate (tmp, param.rot);
    final_pos2 = translate (tmp, param.trans);
    MI = mutual_information (LFM1, final_pos2, LFM2, param, 0);
    return;
end

%% coarse alignment of LFM2 to LFM1

% rotate
param.centroid = calc_centroid(LFM2,param);
canonical = translate (pos2, -param.centroid);
rotated = rotate (canonical,param.angle);
param.rot = param.angle;
if param.just_MI==false
new = translate (rotated, param.centroid);
param.trans = param.centroid;
else
    new = translate (rotated, param.trans);
end

if param.just_MI==false
    if isempty(param.offset)
        % estimate offsets
        offsets = estimate_offsets(LFM1, LFM2, new, param);
%         for i=1:3
%             if ~isempty(param.offset[i])
%                 offsets(i) = param.offset[i];
%             end
%         end
        param.offset = offsets;
    else
        offsets = param.offset;
    end
    fprintf('\nestimated offsets = [%f %f %f]\n',offsets(1),offsets(2),offsets(3));
    new = translate (new, offsets);
    param.trans = param.trans + offsets;
end

%keyboard

%% plot data
if param.plot
    param = save_1d_max_projections(LFM1, LFM2, new, param,...
        '1d_max_projections_presim');
    %drawnow
    %keyboard
    save_2d_max_projections(LFM1, LFM2, new, param, 0,...
        '2d_max_projections_presim');
    %drawnow
    %keyboard
    save_2d_max_projections_compact(LFM1, LFM2, new, param, 0,...
        '2d_max_projections_compact_presim');
    %drawnow
    %keyboard
    save_2d_contour_plots(LFM1, LFM2, pos1, pos2, new, param,...
        'contour_plots_presim');
    drawnow
    %keyboard
end

%% null distribution
% for random rotations and offsets (within some range)
% calculate mutual information, MI
[cdf, centers, nullMIvec, param] = null_distribution (LFM1, LFM2, canonical, param);
param.cdf = cdf;
param.centers = centers;
param.nullMIvec = nullMIvec;

%keyboard
%% print MI of volumes pre- and post-coarse registration
% pre-coarse
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, pos2, LFM2, param, 0);
else
    disp('WTF!');
    keyboard;
end
w = find(param.centers>MI,1); % keep only first instance
if MI > param.bestMI
    fprintf('\npre-coarse MI frac = %.5g\n\n',double(MI)/double(param.bestMI));
elseif ~isempty(w)
    fprintf('\npre-coarse MI frac = %.5g\n\n',param.cdf(w));
else
    disp('WTF?');
    keyboard
end
% post-coarse
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, new, LFM2, param, 0);
    if MI == 0
        disp('WTF?');
        keyboard
    end
else
    disp('WTF!');
    keyboard;
end
w = find(param.centers>MI,1); % keep only first instance
if MI > param.bestMI
    fprintf('\npost-coarse MI frac = %.5g\n\n',double(MI)/double(param.bestMI));
elseif ~isempty(w)
    fprintf('\npost-coarse MI frac = %.5g\n\n',param.cdf(w));
else
    disp('WTF?');
    keyboard
end

if param.just_MI==true
    fprintf('MI==%d\n',MI);
    param
    elapsedTime = toc
    diary off;
    return;
end


%% simulated annealing
[new, param] = simulated_annealing (LFM1, new, LFM2, canonical, param);
if param.plot
    save_stats(param);
end

%keyboard

%% combine registered volumes
if param.savevol
    % would be nice to re-load full unclipped volume here
    % but doing so would jack up the next few plots.
    %comb = combineVols (LFM1, LFM2, param);
    comb = combineVols_iter_dim3 (LFM1, LFM2, param);
end

if param.plot
    param = save_1d_max_projections(LFM1, LFM2, new, param,...
        '1d_max_projections_postsim');
    %drawnow
    %keyboard
    save_2d_max_projections(LFM1, LFM2, comb, param, 1, ...
        '2d_max_projections_postsim');
    %drawnow
    %keyboard
    save_2d_max_projections_compact(LFM1, LFM2, comb, param, 1, ...
        '2d_max_projections_compact_postsim');
    %drawnow
    %keyboard
    save_2d_contour_plots(LFM1, LFM2, pos1, pos2, new, param, ...
        'contour_plots_postsim');
    drawnow
    %keyboard
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
if 0>1
    % assume 16 bit
    dr = 16;
else
    % do NOT assume 16 bit
    m1 = max(max(max(LFM1)));
    m2 = max(max(max(LFM2)));
    dr = ceil(log2(single(max(m1,m2))));
end
%edges = linspace(0,2^dr,round(dr/16*param.N));
edges = [0:2^dr]-0.5;
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

% if i==1
%     param.threshold1 = 0.1;
% else
%     param.threshold1 = edges(i);
% end
% if j==1
%     param.threshold2 = 0.1;
% else
%     param.threshold2 = edges(j);
% end
% param.threshold1 = centers(i);
% param.threshold2 = centers(j);
param.threshold1 = edges(i);
param.threshold2 = edges(j);
fprintf('threshold1 = %f\n',param.threshold1);
fprintf('threshold2 = %f\n',param.threshold2);
param.contour_int1 = param.threshold1;
param.contour_int2 = param.threshold2;
% param.contour_int1 = centers(i);
% param.contour_int2 = centers(j);
fprintf('contour_int1 = %f\n',param.contour_int1);
fprintf('contour_int2 = %f\n\n\n',param.contour_int2);
if param.plot & ~param.justCalcMI
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
    xlim([0 10*centers(i)]);
    mystr1 = [sprintf('%3.3f%% of voxels in LFM1 exceed %.1f in intensity\n',...
        100*(1-cdf_LFM1(i)), param.threshold1 )...
              sprintf('%3.3f%% of voxels in LFM2 exceed %.1f in intensity\n',...
        100*(1-cdf_LFM2(j)), param.threshold2 )...
        sprintf('\nLFM1 has %d voxels\n', numel(LFM1))...
        sprintf('LFM2 has %d voxels', numel(LFM2))...
        ];
   
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1);
    mystr2 = [num2str(round(dr/16*param.N)) ' bins in distribution'];
    text(0.4,0.3,mystr2,'FontSize',9,'Color',[0 0 0],'Interpreter','none');
    text(0.4,0.5,mystr1,'FontSize',9,'Color',[0 0 0],'Interpreter','none');
    
    %keyboard
    
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

% if i == 1
%   % i equaling 1 leads to a bug in later calls
%   % just use the next largest value in the cdf
%   % fixes a problem with large sparse images where < 10^-4
%   % fraction of the pixels are non-zero
%   i = 2
% end

end


function out = combineVols_iter_dim3 (LFM1, LFM2, param)
out = zeros(size(LFM2),'uint16');
% chunk LFM2 in third (z) dimension
% full frame in first (y) and second (x) dimensions
%M = 20; % number of chunks
s = size(LFM2);
%N = ceil(s(3)/M); % number of voxels/frames per chunk
F = s(1)*s(2); % number of voxels per frame
fprintf('\n');
% for each slice of dim 3 in LFM2
if param.parallel
    parfor i = 1:s(3)
        %disp(sprintf('chunk %d of %d',i,M));
        % calculate index range
        linind = [(i-1)*F+1  i*F];
        if linind(1)>numel(LFM2)
            disp('LFM2 is complete!');
            continue;
        end
        if linind(2)>numel(LFM2)
            linind(2) = numel(LFM2);
        end
        fprintf('zindex %2d of %2d, linear index %12d to %12d\n',i,s(3),linind(1),linind(2));
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
        out(:,:,i) = combine (LFM1, LFM2, final_pos2, linind, param);
    end
else
    for i = 1:s(3)
        %disp(sprintf('chunk %d of %d',i,M));
        % calculate index range
        linind = [(i-1)*F+1  i*F];
        if linind(1)>numel(LFM2)
            disp('LFM2 is complete!');
            continue;
        end
        if linind(2)>numel(LFM2)
            linind(2) = numel(LFM2);
        end
        fprintf('zindex %2d of %2d, linear index %12d to %12d\n',i,s(3),linind(1),linind(2));
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
        out(:,:,i) = combine (LFM1, LFM2, final_pos2, linind, param);
    end
end
if param.saveTIF
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
    elseif exist('A','var')
        if isa(A,'uint16')
            XguessSAVE1 = uint16(A);
        elseif isa(A,'double')
            XguessSAVE1 = uint16(A);
        elseif isa(A,'logical')
            %XguessSAVE1 = (2^16-1)*uint16(A);
            XguessSAVE1 = 2^8*uint16(A);
        else
            disp('Warning input data is not 16 bit.');
            keyboard
        end
        clear A;
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
% assume param.interp = 1 or even
if ~(param.interp==1 || mod(param.interp,2)==0)
    fprintf('Error. param.interp is assumed to 1 or an even number (2,4,6,etc)\n.');
    fprintf('param.interp = %d\n',param.interp);
    exit
end
if param.interp>1
    % initialize container of new size
    s = size(LFM);
    out = zeros(s(1),s(2),param.interp*s(3),'uint16');
    boundary = param.interp/2 + 2;
    a = 1;
    b = 0;
    j = 1;
    A = LFM(:,:,j);
    B = LFM(:,:,j+1);
    del = 1/param.interp;
    N = param.interp*s(3)+1;
    last_v = a*A+b*B;
    for i=1:N
        if i < boundary
            if i>1
                out(:,:,i-1) = LFM(:,:,1);
            end
        elseif i > (N-(boundary-1))
            out(:,:,i-1) = LFM(:,:,end);
        else
            v = a*A+b*B;
            out(:,:,i-1) = (last_v + v)/2;
            last_v = v;
            a = a-del;
            b = b+del;
            if a==0
                a=1;
                b=0;
                j=j+1;
                A = LFM(:,:,j);
                B = LFM(:,:,j+1);
            end
        end
    end
else
    out = LFM;
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
s = size(LFM2);
out = reshape(out,s(1),s(2));
%max(max(max(out)))
end


function [cdf, centers, nullMIvec, param] = null_distribution (LFM1, LFM2, canonical, param)
nullf = [param.savePath param.inputFileName{1}(1:end-4) '_null.mat'];
if exist(nullf,'file') == 2
    load(nullf,'nullMIvec','bestMI','bestd','bestr');
    param.bestMI = bestMI;
    param.bestd = bestd;
    param.bestr = bestr;
else
    nullMIvec = [];
    param.bestMI = 0;
    param.bestd = [];
    param.bestr = [];
end
LL = length(nullMIvec);
fprintf('\nnull distribution has N = %d\n',LL);
% determine gain so offset limit = half or quarter of volume limits
% and so that rotation limit = pi
a = max ([ size(LFM1) size(LFM2)]);
offset_limit = 0.15 * a * param.voxel_y;
tmpt = param.trans_amp;
param.trans_amp = offset_limit;
%tmp = param.rot_amp;
gain = pi;
tmpr = param.rot_amp;
if param.confocal
    param.rot_amp = [1 1 1];
else
    param.rot_amp = [1 0 0];
end
%profile on;
fprintf('\nCount    Mutual_Information                Offset [um]                        Rotation [radians]\n');
% if param.parallel
%     parfor i=1:param.Nnull
%         MI = 0;
%         % perturb pos
%         % randomly pick a translation vector and rotation vector
%         % to be added to current location
%         [d,r] = perturb(param,gain);
%         % apply transformation
%         % rotate an amount r PLUS param.rot
%         % thus param.rot tracks the current rotation
%         %rotated = rotate (canonical,r);
%         % translate rotated by an amount param.trans+d
%         % thus param.trans tracks the current position
%         new = translate (rotate (canonical, r), param.centroid+d);
%         %new = translate (rotated, param.centroid+d);
%         % measure mutual_information
%         if strcmp(param.myfunc_MI,'multiply')
%             MI = mutual_information (LFM1, new, LFM2, param);
%             %     elseif strcmp(param.myfunc_MI,'multiply_sqrt')
%             %         MI = mutual_information_sqrt (LFM1, new, LFM2, param);
%         else
%             disp('WTF!');
%             keyboard;
%         end
%         nullMIvec = [nullMIvec MI];
%         str = sprintf('i = %4d, MI = %16.0f, d = [%7.3f0 %7.3f0 %7.3f0], r = [%7.3f0 %7.3f0 %7.3f0]',i,MI,d(1),d(2),d(3),r(1),r(2),r(3));
%         disp(str);
%         %if MI>bestMI
%         %    bestMI = MI;
%         %    bestd = d;
%         %    bestr = r;
%         %end
%         %         profile off
%         %         profile viewer
%     end
% else
    for i=1:param.Nnull
        MI = 0;
        % perturb pos
        % randomly pick a translation vector and rotation vector
        % to be added to current location
        [d,r] = perturb(param, gain, gain);
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
            MI = mutual_information (LFM1, new, LFM2, param, 0);
            %     elseif strcmp(param.myfunc_MI,'multiply_sqrt')
            %         MI = mutual_information_sqrt (LFM1, new, LFM2, param);
        else
            disp('WTF!');
            keyboard;
        end
        nullMIvec = [nullMIvec MI];
        str = sprintf('i = %4d, MI = %16.0f, d = [%7.3f0 %7.3f0 %7.3f0], r = [%7.3f0 %7.3f0 %7.3f0]',i,MI,d(1),d(2),d(3),r(1),r(2),r(3));
        disp(str);
        if MI>param.bestMI
           param.bestMI = MI;
           param.bestd = d;
           param.bestr = r;
        end
        %         profile off
        %         profile viewer
    end
%end
param.rot_amp = tmpr;
param.trans_amp = tmpt;

% add to existing null distribution if any
fprintf('\nnull distribution has N = %d, was N = %d\n',length(nullMIvec),LL);
bestMI = param.bestMI;
bestd = param.bestd;
bestr = param.bestr;
save(nullf,'nullMIvec','bestMI','bestd','bestr');
fprintf('\nnull distribution saved as %s\n',nullf);

if ~param.parallel
    fprintf('\nBest of null:\n');
    fprintf('MI = %d\n',param.bestMI);
    fprintf('d = [%f %f %f]\n',param.bestd(1),param.bestd(2),param.bestd(3));
    fprintf('r = [%f %f %f]\n\n',param.bestr(1),param.bestr(2),param.bestr(3));
end

% CDF
centers = linspace(min(nullMIvec),max(nullMIvec),param.Nnull);
del = 0.5 * (centers(2)-centers(1));
edges  = [centers-del centers(end)+del];
%centers = 0.5 * (edges(1:end-1)+edges(2:end));
% this works. to check: sum(h_LFM1) == numel(LFM1)
h_MI = histcounts(nullMIvec,edges);
% normalize cdf to 1
nh_MI = h_MI / sum(h_MI);

total = 0;
cdf = [];
for i=1:length(nh_MI)
    total = total + nh_MI(i);
    cdf = [cdf total];
end

if param.plot
    % plot PDF
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
    set(gca,'yscale','log');
    
    str=sprintf('%s%s_null_distribution.png',param.savePath,param.timestamp);
    save_plot(f, str);
    
    
    f = figure;
    %plot(centers,cdf);
    %semilogx(centers,cdf);
    if centers(1)>0
        loglog(centers,cdf);
    else
        loglog([1e0 centers(2:end)],cdf);
    end
    if strcmp(param.myfunc_MI,'multiply')
        xlabel('mutual information = sum(LFM1*LFM2)');
    elseif strcmp(param.myfunc_MI,'multiply_sqrt')
        xlabel('mutual information = sum(sqrt(LFM1*LFM2))');
    else
        disp('WTF!');
        keyboard;
    end
    ylabel('normalized cumulative density function');
    %ylim([0 1]);
    title(['null distribution (bootstrapped from'...
        sprintf(' %d random registrations)',length(nullMIvec))]);
    
    str=sprintf('%s%s_null_cdf.png',param.savePath,param.timestamp);
    save_plot(f, str);
end
end

function p = estimateP (param, canonical, LFM1, LFM2, new, T, gain)
init_MI = mutual_information (LFM1, new, LFM2, param, 0);
% Pvec stores history of moves (1) and no moves (0)
Pvec = [];
%N = 10 / param.Pmelt;
j = param.max_moves;
hot = 1;
while j > 0 && numel(Pvec) < param.min_moves
    %fprintf('estimateP: %d of %d\n',j,N);
    % perturb pos
    % randomly pick a translation vector and rotation vector
    % to be added to current location
    [d,r] = perturb(param, gain, gain);
    % hot is an index in range 1 to 6 that selects one dimension out of
    % six (three for translation in dim 1,2,3 and three for rotation about dim 1,2,3)
    if hot<4
        r = zeros(1,3);
        t = r;
        t(hot) = d(hot);
        d = t;
    else
        d = zeros(1,3);
        t = d;
        t(hot-3) = r(hot-3);
        r = t;
    end
    hot = hot+1;
    if hot>6
        hot = 1;
    end
    %str0 = sprintf('d = [%7.3g0 %7.3g0 %7.3g0], r = [%7.3g0 %7.3g0 %7.3g0]',d(1),d(2),d(3),r(1),r(2),r(3));
    % apply transformation
    % rotate an amount r PLUS param.rot
    % thus param.rot tracks the current rotation
    rotated = rotate (canonical,param.rot+r);
    % translate rotated by an amount param.trans+d
    % thus param.trans tracks the current position
    new = translate (rotated, param.trans+d);
    %str1 = print_param(param);
    % measure mutual_information
    if strcmp(param.myfunc_MI,'multiply')
        MI = mutual_information (LFM1, new, LFM2, param, 0);
    else
        disp('WTF!');
        keyboard;
    end
    delmi = MI - init_MI;
    %str2 = sprintf('test MI = %7.3g, delmi = %7.3g',MI,delmi);
    %str3 = 'MI increased';
    % if mutual information dropped
    if delmi < 0.0
        p = exp(delmi/T);
        rnd = rand(1);
        
        if rnd < p
            % accept decrease in MI mutual information
            %str3 = sprintf('Accept %.3g < %.3g',rnd,p);
            Pvec = [Pvec 1];
        else
            % reject move
            Pvec = [Pvec 0];
            %str3 = sprintf('Reject %.3g > %.3g',rnd,p);
        end
%     else
%         disp('delmi=>0');
%         keyboard;
    end
    j=j-1;
    %str4 = sprintf('%d %d, T = %1.0e',i,j,T);
    %fprintf('%7s%75s  %84s  %40s  %22s\n',str4,str0,str1,str2,str3);
end
if numel(Pvec) == 0
    str1 = 'Error. The goal of this function is to measure the frequency';
    str2 = 'at which random moves that decrease mutual information';
    str3 = '(ie worsen registration) are accepted randomly.';
    str4 = sprintf('However, after %d moves, every move increased the MI!',N);
    str5 = 'Basically, this function assumes a good coarse registration';
    str6 = 'so that some moves do worsen the registration';
    fprintf('%s\n%s\n%s\n%s\n%s\n%s\n',str1,str2,str3,str4,str5,str6);
    keyboard
    %exit;
end
fprintf('Of %d moves, %d decreased MI, and of those %d were accepted\n',param.max_moves-j,numel(Pvec),sum(Pvec));

p = sum(Pvec)/numel(Pvec);
end


function T = find_melting_T (LFM1, new, LFM2, canonical, param, gain)
if param.T_fast
    T = param.T0;
    return;
end
p = [];
for t = param.Trange
    p = [p estimateP(param, canonical, LFM1, LFM2, new, t, gain) ];
end

if param.Pmelt < p(1) || param.Pmelt > p(end)
    if param.Pmelt < p(1)
        a = param.Trange;
        fprintf('Warning. Temperature range [%f %f] corresponds to acceptance probabilities [%f %f].\n',a(1),a(end),p(1),p(end));
        fprintf('Temperature range is too high to reach Pmelt = %f.\n',param.Pmelt);
        tvec = [a(1)*a(1)/a(end)];
        while tvec(end)*10 <= a(1)
            tvec = [tvec tvec(end)*10];
        end
        fprintf('Trying a lower temperature range [%f %f].\n',tvec(1),tvec(end));
        p = [];
        for t = tvec
            p = [p estimateP(param, canonical, LFM1, LFM2, new, t, gain) ];
        end
    elseif param.Pmelt > p(end)
        a = param.Trange;
        fprintf('Warning. Temperature range [%f %f] corresponds to acceptance probabilities [%f %f].\n',a(1),a(end),p(1),p(end));
        fprintf('Temperature range is too low to reach Pmelt = %f.\n',param.Pmelt);
        tvec = [a(end)*a(end)/a(1)];
        while tvec(1) >= a(end)
            tvec = [tvec(1)/10 tvec];
        end
        fprintf('Trying a higher temperature range [%f %f].\n',tvec(1),tvec(end));
        p = [];
        for t = tvec
            p = [p estimateP(param, canonical, LFM1, LFM2, new, t, gain) ];
        end
    end
    if param.Pmelt < p(1) || param.Pmelt > p(end)
        fprintf('Error. Temperature range [%f %f] corresponds to acceptance probabilities [%f %f].\n',a(1),a(end),p(1),p(end));
        fprintf('Temperature range is still incorrect to reach Pmelt = %f.\n',param.Pmelt);
        exit
    end
end

Fit = polyfit(log10(param.Trange),p,1);
log10T = (param.Pmelt-Fit(2))/Fit(1);
T = 10^log10T;
end

function [new, param] = simulated_annealing (LFM1, new, LFM2, canonical, param)
param.transvec = [param.trans];
param.offsetvec = [param.offset];
param.rotvec = [param.rot];
param.Pvec = [];
%param.trans_amp = 3.5; %um, half diameter of neuron
% calc rotational gain
a = size(LFM2);
half_span = 0.5*[a(1)*param.voxel_y a(2)*param.voxel_x];
radius = sqrt( sum( half_span .* half_span ) );
gain = param.trans_amp / radius;
gain0 = gain;
% calculate MI
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (LFM1, new, LFM2, param, 0);
else
    disp('WTF!');
    keyboard;
end
%last_MI = MI;
param.MIvec = MI;
% CDF
param.cdfvec = [];
w = find(param.centers>param.MIvec(end),1); % keep only first instance
if ~isempty(w)
    param.cdfvec = [param.cdfvec param.cdf(w)];
else
    val = double(param.MIvec(end))/double(param.bestMI);
    param.cdfvec = [param.cdfvec val] ;
end
% set initial T
% Start with T sufficiently high to "melt" the system
% later to guarantee melting, if needed T will be increased until
% P (accepting a decrease in MI) = 60%
T = find_melting_T(LFM1, new, LFM2, canonical, param, gain);
param.Tvec = T;
%param = set_prate (MI,param);
%p = param.init_p;
% set max number of temperature changes and mean changes
Tchanges = param.TC0;
if strcmp(param.anneal,'exp')
    param.Trate = 10^(log10(1e-3)/Tchanges);
elseif  strcmp(param.anneal,'linear')
    param.Trate = T / Tchanges;
else
    disp('WTF!');
    keyboard;
end
fprintf('\n\nTrate = %f\n\n',param.Trate);

param.Dvec = [];
% while system not frozen and more temperature changes are allowed
% profile on;
disp('Count      Perturbation     Transformation      Overlap     Probability,Decision   ');
%%%% TODO -Add condition 'if no change in voxel overlap between LFM1 + DDLFM
i = 1;
pass = 0;
hot = 1;
hotvec = [1:6];
hotvec = hotvec(randperm(numel(hotvec)));
phase = 2;
deleteme = false;
last_pass = 0;
while 1 > 0
    if phase==2 && (pass>param.frozen2 || i>param.frozen3 || gain < gain / param.gain_limit)
        fprintf('(\n\n%d passes in a row. The model is frozen.\n\n',param.frozen2);
        break;
    end
    if phase==1 && pass>param.frozen1
        fprintf('(\n\n%d passes in a row. Switching to fine search.\n\n',param.frozen1);
        phase = 2;
        pass = 0;
    end
    i = i+1;
    if strcmp(param.anneal,'exp')
        T = T * param.Trate;
    elseif  strcmp(param.anneal,'linear')
        T = T - param.Trate;
    else
        disp('WTF!');
        keyboard;
    end
%     if i > Tchanges
%         param.Tvec = [param.Tvec T];
%         break;
%     end
    if T < param.Tmin
        T = param.Tmin;
%         param.Tvec = [param.Tvec T];
%         break;
    end
    if pass > (last_pass + param.last_pass_N)
        gain = gain / param.gain_scale;
        last_pass = pass;
        str9 = sprintf('gain reduced from %f to %f', gain * param.gain_scale, gain);
    else
        str9 = '';
    end
    % perturb pos
    % randomly pick a translation vector and rotation vector
    % to be added to current location
    [d,r] = perturb(param,gain, gain0);
    if phase == 2
        rhot = hotvec(hot);
        if rhot<4
            r = zeros(1,3);
            t = r;
            t(rhot) = d(rhot);
            d = t;
        else
            d = zeros(1,3);
            t = d;
            t(rhot-3) = r(rhot-3);
            r = t;
        end
        hot = hot+1;
        if hot>6
            hot = 1;
            hotvec = hotvec(randperm(numel(hotvec)));
        end
    end
    str0 = sprintf('d = [%7.3f %7.3f %7.3f], r = [%7.5f %7.5f %7.5f]',d(1),d(2),d(3),r(1),r(2),r(3));
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
        MI = mutual_information (LFM1, new, LFM2, param, 0);
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
        pass = 0;
        last_pass = 0;
    elseif delmi == 0
        str3 = 'Reject - MI unchanged';
        pass = pass + 1;
        param.rot = param.rot - r;
        param.trans = param.trans - d;
        tmp = param.MIvec(end);
        param.MIvec = [param.MIvec tmp];
        %           str3 = 'Accept - MI unchanged';
        %           pass = pass + 1;
        %         MI2 = mutual_information (LFM1, new, LFM2, param, 1);
        %         param.rot = param.rot - r;
        %         param.trans = param.trans - d;
        %         rotated1 = rotate (canonical,param.rot);
        %         new1 = translate (rotated1, param.trans);
        %         MI1 = mutual_information (LFM1, new1, LFM2, param,1);
        %         keyboard
        
    else
        %         # else accept D with probability P = exp(-E/kBT)
        %         # using (psuedo-)random number uniformly distributed in the interval (0,1)
        p = exp(delmi/T);
        param.Pvec = [param.Pvec p];
        % Specifically, if random number is less the P, then accept
        rnd = rand(1);
        if rnd < p
            % accept decrease in MI mutual information
            str3 = sprintf('Accept %.3g < %.3g',rnd,p);
            param.MIvec = [param.MIvec MI];
            pass = 0;
            last_pass = 0;
        else
            % reject move
            str3 = sprintf('Reject %.3g > %.3g',rnd,p);
            param.rot = param.rot - r;
            param.trans = param.trans - d;
            tmp = param.MIvec(end);
            param.MIvec = [param.MIvec tmp];
            pass = pass + 1;
        end
    end
    str4 = sprintf('%d, pass = %d, T = %g',i,pass, T);
    str6 = sprintf('Final MI = %d',param.MIvec(end));
    str5 = print_param(param);
    val7 = param.trans - param.centroid;
    str7 = sprintf('final offset = [%6.6g0 %6.6g0 %6.6g0]',val7(1),val7(2),val7(3));
    w = find(param.centers>param.MIvec(end),1); % keep only first instance
    if param.MIvec(end) >= param.bestMI
        val = double(param.MIvec(end))/double(param.bestMI);
        str8 = sprintf('MI frac = %.5g, siman exceeds null',val);
        param.cdfvec = [param.cdfvec val] ;
    elseif ~isempty(w)
        str8 = sprintf('MI frac = %.5g, null exceeds siman',param.cdf(w));
        param.cdfvec = [param.cdfvec param.cdf(w)];
    else
        disp('WTF?');
        keyboard
    end
    fprintf('%7s%75s  %84s  %40s  %22s  %84s  %20s %40s %20s %s\n',str4,str0,str1,str2,str3,str5,str6,str7,str8,str9);
    param.Tvec = [param.Tvec T];
    param.transvec = [param.transvec; param.trans];
    param.offsetvec = [param.offsetvec; val7];
    param.rotvec = [param.rotvec; param.rot];
    param.Dvec = [param.Dvec; d];
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


% plot Tvec
h = figure;
plot(1:numel(param.Tvec),param.Tvec);
xlabel('iteration');
ylabel('Temperature (-)');
title('simulated annealing schedule');
if 0>1
    fname = sprintf('%s_p.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_T.png',prefix);
    save_plot(h, str);
end

% plot MIvec
h = figure;
plot(1:numel(param.MIvec),log10(param.MIvec));
hold on;
plot(xlim,log10(param.bestMI)*ones(2,1),'--r');
hold off;
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

% plot CDFvec
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

% plot amplitude of 3D random walk
h = figure;
a = size(param.Dvec);
subplot(1,3,1);
plot(1:a(1),param.Dvec(:,1));
xlabel('iteration');
ylabel('random walk in dim 1 (um)');
hold on;
plot(1,param.Dvec(1,1),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,2);
plot(1:a(1),param.Dvec(:,2));
xlabel('iteration');
ylabel('random walk in dim 2 (um)');
hold on;
plot(1,param.Dvec(1,2),'ro');
hold off;
xlim([1 a(1)]);

subplot(1,3,3);
plot(1:a(1),param.Dvec(:,3));
xlabel('iteration');
ylabel('random walk in dim 3 (um)');
hold on;
plot(1,param.Dvec(1,3),'ro');
hold off;
xlim([1 a(1)]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'Attempted translations';
text(0.4,0.97,str,'FontSize',12,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_d.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_d.png',prefix);
    save_plot(h, str);
end

% plot offset
h = figure;
a = size(param.offsetvec);
subplot(1,3,1);
plot(1:a(1),param.offsetvec(:,1));
xlabel('iteration');
ylabel('offset in dim 1 (um)');
hold on;
plot(1,param.offsetvec(1,1),'ro');
plot(a(1),param.bestd(1),'gs');
hold off;
xlim([1 a(1)]);

subplot(1,3,2);
plot(1:a(1),param.offsetvec(:,2));
xlabel('iteration');
ylabel('offset in dim 2 (um)');
hold on;
plot(1,param.offsetvec(1,2),'ro');
plot(a(1),param.bestd(2),'gs');
hold off;
xlim([1 a(1)]);

subplot(1,3,3);
plot(1:a(1),param.offsetvec(:,3));
xlabel('iteration');
ylabel('offset in dim 3 (um)');
hold on;
plot(1,param.offsetvec(1,3),'ro');
plot(a(1),param.bestd(3),'gs');
hold off;
xlim([1 a(1)]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
str = 'trajectory of simulated annealing (Green squares are best of null)';
text(0.4,0.97,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_offset.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_offset.png',prefix);
    save_plot(h, str);
end

% plot rot
h = figure;
plot_rot(h,1,param);
plot_rot(h,2,param);
plot_rot(h,3,param);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
%str = 'trajectory of simulated annealing (RHS axis is degrees, green squares are best of null.)';
str = 'trajectory of simulated annealing (RHS axis is degrees)';
text(0.1,0.98,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');

if 0>1
    fname = sprintf('%s_rot.fig',prefix);
    savefig(h,fname);
else
    str=sprintf('%s_rot.png',prefix);
    save_plot(h, str);
end

end

function plot_rot (h, i, param)
figure(h);

a = size(param.rotvec);
x = 1:a(1);
yrad = param.rotvec(:,i);

subplot(1,3,i);
yyaxis left;
plot(x,yrad),'k';
hold on;
plot(1,param.rotvec(1,i),'ro');
%plot(a(1),param.bestr(i),'gs');
hold off;

xlabel('iteration');
ylabel(['rotation around dim ' num2str(i) ' (radians)']);
xlim([1 a(1)]);
y = ylim;

scale = 180 / pi;
%ydeg = yrad * scale;
yyaxis right;
%h2 = plot(x,ydeg,'k');
ny = y*scale;
ylim(ny);
%set(h2, 'Visible' ,'off');
end




function out = print_param (param)
a=param.trans;
b=param.rot;
out = sprintf('trans = [%6.6f %6.6f %6.6f], rot = [%7.6f %7.6f %7.6f]',a(1),a(2),a(3),b(1),b(2),b(3));
end


function mi = mutual_information (LFM1, new, LFM2, param, deleteme)
scale = [1/param.voxel_y 1/param.voxel_x 1/param.voxel_z];
new_pixels = ceil(new.*scale);
s = size(LFM1);
% find rows of new_pixels that overlap LFM1 and save as index
index = find( new_pixels(:,1)<=s(1) & new_pixels(:,2)<=s(2) & new_pixels(:,3)<=s(3) ...
    & new_pixels(:,1)>0 & new_pixels(:,2)>0 & new_pixels(:,3)>0 );
% in LFM1, lookup intensity value at 
i1 = LFM1(sub2ind(s,new_pixels(index,1),new_pixels(index,2),new_pixels(index,3)));
%i1 = LFM1(new_pixels(index,1),new_pixels(index,2),new_pixels(index,3));
i2 = LFM2(param.index2(index));
mi = sum(double(i1).*double(i2));
if deleteme
    keyboard
    %save('before.mat','new','index','i1','i2','mi');
end
end


function out = rotation_matrix (angle)
a = [1 0 0;...
    0 cos(angle(1)) -sin(angle(1));...
    0 sin(angle(1)) cos(angle(1))];
b = [cos(angle(2)) 0 sin(angle(2));...
    0 1 0;...
    -sin(angle(2)) 0 cos(angle(2))];
c = [cos(angle(3)) -sin(angle(3)) 0;...
    sin(angle(3)) cos(angle(3)) 0;...
    0 0 1];
%out = a*b*c;
out = c*b*a;
end


function [out] = rotate (pos, angle)
rot = single( rotation_matrix (angle) );
%out = pos*rot;
out = rot*pos';
out = out';
end


function [out] = translate (pos, delta)
D = ones(size(pos),'single').*delta;
out = pos + D;
end

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


function [d,r] = perturb (param, gain, gain0)
% but gain is not used for d
% so introduce gain ratio to turn down d
a = gain / gain0;
d1 = random('unif',-1,1) * param.trans_amp * a;
d2 = random('unif',-1,1) * param.trans_amp * a;
d3 = random('unif',-1,1) * param.trans_amp * a;
% how gain was originally used
r1 = random('unif',-1,1) * param.rot_amp(1) * gain;
r2 = random('unif',-1,1) * param.rot_amp(2) * gain;
r3 = random('unif',-1,1) * param.rot_amp(3) * gain;
d = [d1 d2 d3];
r = [r1 r2 r3];
%r = [r1 0 0];
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
[~, index] = min(abs(u_index - numel(u_acor)/2));
u_index = u_index(index);
v_index = find(v_acor == max(v_acor));
[~, index] = min(abs(v_index - numel(v_acor)/2));
v_index = v_index(index);
w_index = find(w_acor == max(w_acor));
[~, index] = min(abs(w_index - numel(w_acor)/2));
w_index = w_index(index);

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



function param = save_1d_max_projections (LFM1, LFM2, new, param, str)
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

if 0 < 1
    T = '      (scaled by total intensity in sample)';
    scale_d1 = [ sum(LFM1_d1) sum(LFM2_d1) sum(new_d1)];
    scale_d2 = [ sum(LFM1_d2) sum(LFM2_d2) sum(new_d2)];
    scale_d3 = [ sum(LFM1_d3) sum(LFM2_d3) sum(new_d3)];
    yl = 'normalized intensity';
    %disp(T);
elseif 0 < 1
    T = '      (scaled by max intensity in each projected volume)';
    scale_d1 = single([ max(LFM1_d1) max(LFM2_d1) max(new_d1)]);
    scale_d2 = single([ max(LFM1_d2) max(LFM2_d2) max(new_d2)]);
    scale_d3 = single([ max(LFM1_d3) max(LFM2_d3) max(new_d3)]);
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
semilogy(LFM1_d1/scale_d1(1),'.','Color',colors{1});
hold on;
semilogy(LFM2_d1/scale_d1(2),'.','Color',colors{2});
semilogy(new_d1/scale_d1(3),'-','Color',colors{4});
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
semilogy(LFM1_d2/scale_d2(1),'.','Color',colors{1});
hold on;
semilogy(LFM2_d2/scale_d2(2),'.','Color',colors{2});
semilogy(new_d2/scale_d2(3),'-','Color',colors{4});
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
semilogy(LFM1_d3/scale_d3(1),'.','Color',colors{1});
hold on;
semilogy(LFM2_d3/scale_d3(2),'.','Color',colors{2});
semilogy(new_d3/scale_d3(3),'-','Color',colors{4});
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

msg = sprintf('LFM1: contour at intensity level %4.2f',param.contour_int1);
text(0.2,0.4,msg,'FontSize',12,'Color',colors{3},'Interpreter','none');

msg = sprintf('LFM2: contour at intensity level %4.2f',param.contour_int2);
text(0.2,0.35,msg,'FontSize',12,'Color',colors{2},'Interpreter','none');

msg = sprintf('LFM2 coarse reg: contour at intensity level %4.2f',param.contour_int2);
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
    if dr == 0 || max(max(max(new)))==0
        fprintf('Warning! dr =0.\n');
        dr = 8;
    end
    imagesc(yz,[0 2^dr]);
    xlabel('three [pixels]');
    ylabel('one [pixels]');
    title('DDLFM');
    %colorbar();
    daspect([1,1,1]);
    set(gca,'XDir','reverse');
    
    xy = squeeze(max(new,[],3));
    subplot(2,6,6);
    dr = ceil(log2(single(max(max(xy)))));
    if dr == 0 || max(max(max(new)))==0
        fprintf('Warning! dr =0.\n');
        dr = 8;
    end
    imagesc(xy,[0 2^dr]);
    xlabel('two [pixels]');
    ylabel('one [pixels]');
    title('DDLFM');
    %colorbar();
    daspect([1,1,1]);
    
    xz = squeeze(max(new,[],1))';
    subplot(2,6,12);
    dr = ceil(log2(single(max(max(xz)))));
    if dr == 0 || max(max(max(new)))==0
        fprintf('Warning! dr =0.\n');
        dr = 8;
    end
    imagesc(xz,[0 2^dr]);
    xlabel('two [pixels]');
    ylabel('three [pixels]');
    title('DDLFM');
    colorbar();
    daspect([1,1,1]);
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.1,0.98,['LFM1 = ' param.inputFilePath1 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
text(0.1,0.96,['LFM2 = ' param.inputFilePath2 param.inputFileName{1}],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
if flag
    text(0.1,0.94,['DDLFM = ' param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
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
    title('DDLFM');
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
    text(0.1,0.94,['DDLFM = ' param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
end
text(0.4,0.1,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
save_plot(f, str);

end
