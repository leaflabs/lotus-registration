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
LFM1 = loadData(f, param, 'LFM1');
f = [param.inputFilePath2 param.inputFileName{1}];
fprintf('\nLoading volume %s\n\n',f);
LFM2 = loadData(f, param, 'LFM2');
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
        fprintf('\nEstimating offsets...\n');
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
        fprintf('\nUsing stored offsets...\n');
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


% The following functions were split into separate files
%   calculate_thresholds
%   derive_threshold
%   combineVols_iter_dim3
%   combineVols
%   loadData
%   interpolate
%   combine
%   null_distribution
%   estimateP
%   find_melting_T
%   simulated_annealing
%   save_plot
%   save_stats
%   plot_rot
%   print_param
%   mutual_information
%   rotation_matrix
%   rotate
%   translate
%   calc_centroid
%   init_pos
%   save_vol
%   perturb
%   print_fraction
%   estimate_offsets
%   save_1d_max_projections
%   save_2d_contour_plots
%   save_2d_max_projections
%   save_2d_max_projections_compact


% Discarded lines were placed into a file named "__discard__.m"