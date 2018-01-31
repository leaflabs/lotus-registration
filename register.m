%parpool(4);
tic

%% load parameters
param = struct;
param = reeset(param);
disp('Loading parameters:');
a = who;
for i=1:length(a)
    mystr = sprintf('param.%s = %s;',a{i},a{i});
    disp(mystr);
    eval(mystr)
end
param.trans_amp = param.scale_trans * param.voxel_x; % um
param.rot_amp = param.scale_rot * pi/800; % radians
param.timestamp = [param.inputFileName{1}(1:end-4) '__' datestr(datetime('now'),'yyyymmdd_HHMMSS')];
fname = sprintf('%s%s.log',param.savePath, param.timestamp);
diary(fname)


%% load volumes
disp(sprintf('\nLoading volume %s\n',param.inputFileName{1}));
f = [param.inputFilePath1 param.inputFileName{1}];
side1 = loadData(f, param);
f = [param.inputFilePath2 param.inputFileName{1}];
side2 = loadData(f, param);
param.voxel_z = param.voxel_z / param.interp;

%% if registration parameters are already known
%  then combine registered volumes
if param.rapid
    % must specify: param.centroid, param.trans, param.rot
    combineVols (side1, side2, param);
    return;
end


%% save TIF versions of input volumes
if param.lfdisplay
    f = param.inputFileName{1};
    outFile = sprintf('%s.tif',f(1:end-4));
    save_vol( side1, param.inputFilePath1, outFile);
    save_vol( side2, param.inputFilePath2, outFile);
end


%% clip volumes for registration
a = size(side1);
if ~isempty(find(param.clip>0))
    disp('Clipping pixels from periphery:');
    disp(sprintf('%d ',param.clip));
    disp(sprintf('\n'));
    side1 = side1( 1+param.clip(1):a(1)-param.clip(2),...
        1+param.clip(3):a(2)-param.clip(4),...
        1+param.clip(5):a(3)-param.clip(6));
end
a = size(side2);
if ~isempty(find(param.clip>0))
    side2 = side2( 1+param.clip(1):a(1)-param.clip(2),...
        1+param.clip(3):a(2)-param.clip(4),...
        1+param.clip(5):a(3)-param.clip(6));
end

%% calculate thresholds
param.N = 10000;
param.xlim_thresh = 0.999; %HARDCODED
param.contour_thresh = 0.99;
param.pop_thresh = 0.95; %HARDCODED
param = calculate_thresholds (side1, side2, param,'cdf_voxel_intensities');

%% voxel positions
param.index1 = find(side1>param.threshold1);
print_fraction(param.index1,side1,'side1');
pos1 = init_pos(param.index1,side1,param);
param.index2 = find(side2>param.threshold2);
print_fraction(param.index2,side2,'side2');
param.size = size(side2);
pos2 = init_pos(param.index2,side2,param);


%% coarse alignment of side2 to side1
param.centroid = calc_centroid(side2,param);
[canonical,param] = translate (pos2, -param.centroid, param);
[rotated,param] = rotate (canonical,param.rot+param.angle, param);
param.rot = param.rot + param.angle;
[new,param] = translate (rotated, param.centroid+param.offset, param);
param.trans = param.centroid + param.offset;


%% plot data
if param.plot
    save_1d_max_projections(side1, side2, pos1, pos2, new, param,...
        '1d_max_projections_presim');
    save_2d_max_projections(side1, side2, new, param, 0,...
        '2d_max_projections_presim');
    save_2d_max_projections_compact(side1, side2, new, param, 0,...
        '2d_max_projections_compact_presim');
    save_2d_contour_plots(side1, side2, pos1, pos2, new, param,...
        'contour_plots_presim');
    drawnow
end

%% null distribution
% for random rotations and offsets (within some range)
% calculate mutual information, MI
param.N = 2000;
[cdf,centers,nullMIvec] = null_distribution (pos1, param.index1, side1, new, side2, canonical, param);
param.cdf = cdf;
param.centers = centers;
param.nullMIvec = nullMIvec;

%% simulated annealing
[new, param] = simulated_annealing (pos1, param.index1, side1, new, side2, canonical, param);
if param.plot
    %plot_pos(pos1, pos2, side1, side2, new, param, 'postsim');
    save_stats(param);
end


%% combine registered volumes
if param.savevol
    % would be nice to re-load full unclipped volume here
    % but doing so would jack up the next few plots.
    comb = combineVols (side1, side2, param);
end

if param.plot
    save_1d_max_projections(side1, side2, pos1, pos2, new, param,...
        '1d_max_projections_postsim');
    save_2d_max_projections(side1, side2, comb, param, 1, ...
        '2d_max_projections_postsim');
    save_2d_max_projections_compact(side1, side2, comb, param, 1, ...
        '2d_max_projections_compact_postsim');
    save_2d_contour_plots(side1, side2, pos1, pos2, new, param, ...
        'contour_plots_postsim');
    drawnow
end

%% save final parameters
param
fname = sprintf('%s%s_parameters.mat',param.savePath,param.timestamp);
save(fname,'param');
elapsedTime = toc
diary off;


%% functions

function param = calculate_thresholds (side1, side2, param, str)
% assume 16 bit
edges = linspace(0,2^16,param.N);
centers = (edges(1:end-1)+edges(2:end))/2;
% this works. to check: sum(h_side1) == numel(side1)
h_side1 = histcounts(side1,edges);
h_side2 = histcounts(side2,edges);

% normalize cdf to 1
h_side1 = h_side1 / sum(sum(sum(h_side1)));
h_side2 = h_side2 / sum(sum(sum(h_side2)));

total = 0;
cdf_side1 = [];
for i=1:length(h_side1)
    total = total + h_side1(i);
    cdf_side1 = [cdf_side1 total];
end
total = 0;
cdf_side2 = [];
for i=1:length(h_side2)
    total = total + h_side2(i);
    cdf_side2 = [cdf_side2 total];
end

f = figure;

plot(centers,cdf_side1);
hold on;
plot(centers,cdf_side2);
xlabel('Intensity of voxel [uint16]');
ylabel('normalized cumulative density function');
title([num2str(param.N) ' bins in distribution']);
legend('side1','side2');
hold off;
i = find(cdf_side1>param.xlim_thresh,1); % keep only first instance
xlim([0 centers(i)]);

i = find(cdf_side1>param.pop_thresh,1); % keep only first instance
j = find(cdf_side2>param.pop_thresh,1); % keep only first instance
p = find(cdf_side1>param.contour_thresh,1); % keep only first instance
q = find(cdf_side2>param.contour_thresh,1); % keep only first instance
mystr = [sprintf('brightest %2.0f%% of voxel population\n',100*(1-param.pop_thresh))...
    'has an intensity larger than'...
    sprintf('\n%3.0f for side1\n',centers(i))...
    sprintf('%3.0f for side2\n\n\n',centers(j))...
    sprintf('side 1 has %d voxels\n',numel(side1))...
    sprintf('side 2 has %d voxels\n',numel(side2))...
    ];
a = xlim;
b=ylim;
text(a(2)/3,b(2)/2,mystr,'FontSize',8,'Color',[0 0 0]);

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);

param.threshold1 = centers(i);
param.threshold2 = centers(j);
param.contour_int1 = centers(p);
param.contour_int2 = centers(q);
end


function out = combineVols (side1, side2, param)
out = zeros(size(side2),'uint16');
% chunk side2 in third (z) dimension
% full frame in first (y) and second (x) dimensions
M = 20; % number of chunks
s = size(side2);
N = ceil(s(3)/M); % number of voxels/frames per chunk
F = s(1)*s(2); % number of voxels per frame
% for each chunk in side2
for i = 1:M
    %disp(sprintf('chunk %d of %d',i,M));
    % calculate index range
    zind = [(i-1)*N+1  i*N];
    if zind(2)>s(3)
        zind(2) = s(3);
    end
    linind = [(i-1)*N*F+1  i*N*F];
    if linind(1)>numel(side2)
        disp('Side2 is complete!');
        break;
    end
    if linind(2)>numel(side2)
        linind(2) = numel(side2);
    end
    disp(sprintf('\nchunk %d of %d, zind %d %d, linind %d %d',i,M, zind(1),zind(2),linind(1),linind(2)));
    % for each voxel in chunk in side2, extract position to final_pos2
    disp('Init_pos');
    final_pos2 = init_pos([linind(1):linind(2)]', side2, param);
    % for each position in final_pos2, apply rot and translation
    disp('Translate');
    [tmp,param] = translate (final_pos2, -param.centroid, param);
    disp('Rotate');
    [tmp,param] = rotate (tmp, param.rot, param); % CHECK
    disp('Translate');
    [final_pos2,param] = translate (tmp, param.trans, param); 
    % for each voxel in side1, find positions in final_pos2 and combine and normalize
    disp('combine');
    out(linind(1):linind(2)) = combine (side1, side2, final_pos2, zind, linind, param);
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
disp(sprintf('Saving combined volume to %s.',outFile));
save(outFile,'XguessSAVE1','-v7.3');
end

function out = loadData (f, param)
load(f);
if exist('Xvolume','var')
    if strcmp(class(Xvolume),'uint16')
        XguessSAVE1 = Xvolume;
    else
        disp('Warning input data is not 16 bit.');
        keyboard
    end
    clear Xvolume;
elseif exist('XguessSAVE1','var')
    disp('XguessSAVE1 found.');
    if ~strcmp(class(XguessSAVE1),'uint16')
        disp('Warning input data is not 16 bit.');
        keyboard
        %disp('Data will be analyzed as 8 bit.');
        %XguessSAVE1 = cast(XguessSAVE1, 'uint16');
    end
else
    disp('WTF?! Unknown data name.');
    return;
end
out = interpolate (XguessSAVE1,param);
end

function out = interpolate (side, param)
% initialize container of new size
s = size(side);
out = zeros(s(1),s(2),param.interp*s(3),'uint16');
boundary = param.interp/2 + 1;
for i=1:param.interp*s(3)
    if i < boundary
        out(:,:,i) = side(:,:,1);
    elseif i > (param.interp*s(3)-(boundary-1))
        out(:,:,i) = side(:,:,end);
    else
        a = round((i-1)/param.interp);
        b = a+1;
        N = 2*param.interp;
        fb = 1+2*mod(i-(param.interp-1),param.interp);
        fa = N-fb;
        out(:,:,i) = fa/N*side(:,:,a) + fb/N*side(:,:,b);
    end
end
end


function out = combine (side1, side2, chunk, zind, linind, param)
% chunk = x,y,z centroids of voxels in side2
% for each centroid in chunk, calculate index of corresponding voxels in side1
L = length(chunk);
scale = [1/param.voxel_y*ones(L,1) ...
    1/param.voxel_x*ones(L,1) ...
    1/param.voxel_z*ones(L,1)];
abc_side1 = ceil(chunk.*scale);
% purge any indices that lie outside of side1
ind = find(abc_side1<1);
[r1,c] = ind2sub(size(abc_side1),ind);
s = size(side1);
ind = find(abc_side1(:,1)>s(1));
[r2,c] = ind2sub(size(abc_side1),ind);
ind = find(abc_side1(:,2)>s(2));
[r3,c] = ind2sub(size(abc_side1),ind);
ind = find(abc_side1(:,3)>s(3));
[r4,c] = ind2sub(size(abc_side1),ind);
rem = unique([r1;r2;r3;r4]);
disp(sprintf('abc_side1 %d, rem %d, r1 %d, r2 %d, r3 %d, r4 %d',...
    length(abc_side1),length(rem),length(r1),length(r2),length(r3),length(r4)));
% chunk is a systematic sweep in z dimensino of side2
% in contrast, abc_side1 are the voxels in side1 coordinate system of a
% rotated and translated side2, so not well behaved
% thus, compute intensity in side2 coordinate system
out = ones(length(chunk),1,'uint16');
%out = zeros(length(chunk),1,'uint16');
%
% how to handle rem entries
% in abc, replace with 1,1,1
abc_side1(rem,:) = [ones(length(rem),1) ones(length(rem),1) ones(length(rem),1)];
%
% in intensity calculation, set all rem entries to zero
ind = sub2ind(size(side1),abc_side1(:,1),abc_side1(:,2),abc_side1(:,3));
i1 = side1(ind);
i1(rem) = 0.0;
i2 = side2(linind(1):linind(2))';

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


function [cdf,centers,nullMIvec] = null_distribution (pos1, index1, side1, new, side2, canonical, param)
nullMIvec = [];
gain = 5;
tmp = param.rot_amp;
param.rot_amp = param.rot_amp * 20;
%profile on;
disp(sprintf('\nCount    Offset                            Rotation                          Mutual_Information'));
for i=1:param.N
    % perturb pos
    % randomly pick a translation vector and rotation vector
    % to be added to current location
    [d,r] = perturb(param,gain);
    str0 = sprintf('i = %d, d = [%7.3f0 %7.3f0 %7.3f0], r = [%7.3f0 %7.3f0 %7.3f0]',i,d(1),d(2),d(3),r(1),r(2),r(3));
    % apply transformation
    % rotate an amount r PLUS param.rot
    % thus param.rot tracks the current rotation
    [rotated,param] = rotate (canonical,r, param);
    % translate rotated by an amount param.trans+d
    % thus param.trans tracks the current position
    [new,param] = translate (rotated, param.centroid+d, param);
    str1 = print_param(param);
    % measure mutual_information
    if strcmp(param.myfunc_MI,'multiply')
        MI = mutual_information (pos1, index1, side1, new, side2, param, 0);
    elseif strcmp(param.myfunc_MI,'multiply_sqrt')
        MI = mutual_information_sqrt (pos1, index1, side1, new, side2, param, 0);
    else
        disp('WTF!');
        keyboard;
    end
    nullMIvec = [nullMIvec MI];
    disp(sprintf('%75s, MI = %f',str0,MI));
    %         profile off
    %         profile viewer
end
param.rot_amp = tmp;
%profile viewer;

% CDF
edges = linspace(0,max(nullMIvec),param.N);
centers = (edges(1:end-1)+edges(2:end))/2;
% this works. to check: sum(h_side1) == numel(side1)
h_MI = histcounts(nullMIvec,edges);
% normalize cdf to 1
h_MI = h_MI / sum(h_MI);

total = 0;
cdf = [];
for i=1:length(h_MI)
    total = total + h_MI(i);
    cdf = [cdf total];
end

% add to existing null distribution if any
nullf = [param.savePath 'null.mat'];
if exist(nullf,'file') == 2
    tmp = nullMIvec;
    load(nullf,'nullMIvec');
    disp(sprintf('nullMIvec was N = %d',length(nullMIvec)));
    nullMIvec = [nullMIvec tmp];
    disp(sprintf('nullMIvec is now N = %d',length(nullMIvec))); 
end
save(nullf,'nullMIvec');

% plot PDF and CDF
f = figure;
histogram(nullMIvec);
if strcmp(param.myfunc_MI,'multiply')
    xlabel('mutual information = sum(side1*side2)');
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    xlabel('mutual information = sum(sqrt(side1*side2))');
else
    disp('WTF!');
    keyboard;
end
ylabel('count');
title(['null distribution (bootstrapped from'...
    sprintf(' %d random registrations)',length(nullMIvec))]);

str=sprintf('%s%s_null_distribution.png',param.savePath,param.timestamp);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);


f = figure;
plot(centers,cdf);
if strcmp(param.myfunc_MI,'multiply')
    xlabel('mutual information = sum(side1*side2)');
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    xlabel('mutual information = sum(sqrt(side1*side2))');
else
    disp('WTF!');
    keyboard;
end
ylabel('normalized cumulative density function');
ylim([0 1]);
title(['null distribution (bootstrapped from'...
    sprintf(' %d random registrations)',length(nullMIvec))]);

str=sprintf('%s%s_null_cdf.png',param.savePath,param.timestamp);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);
end




function [new, param] = simulated_annealing (pos1, index1, side1, new, side2, canonical, param)
%print_param(param);
if strcmp(param.myfunc_MI,'multiply')
    MI = mutual_information (pos1, index1, side1, new, side2, param, 0);
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    MI = mutual_information_sqrt (pos1, index1, side1, new, side2, param, 0);
else
    disp('WTF!');
    keyboard;
end
last_MI = MI;
param = setT0 (MI,param);
param.MIvec = [MI];
% for plotting only
if strcmp(param.myfunc_MI,'multiply')
    MIt = mutual_information (pos1, index1, side1, new, side2, param, param.threshold_plot);
elseif strcmp(param.myfunc_MI,'multiply_sqrt')
    MIt = mutual_information_sqrt (pos1, index1, side1, new, side2, param, param.threshold_plot);
else
    disp('WTF!');
    keyboard;
end
param.MItvec = [MIt];
% set initial T. Start with T sufficiently high to "melt" the system
T = param.T0;
p = param.init_p;
param.Pvec = [p];
% set max number of temperature changes and mean changes
Tchanges = param.TC0;
% while system not frozen and more temperature changes are allowed
%profile on;
param.transvec = [param.trans];
param.rotvec = [param.rot];
% CDF
param.cdfvec = [];
w = find(param.centers>last_MI,1); % keep only first instance
if ~isempty(w)
    param.cdfvec = [param.cdfvec param.cdf(w)];
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
        [rotated,param] = rotate (canonical,param.rot+r, param);
        param.rot = param.rot + r;
        % translate rotated by an amount param.trans+d
        % thus param.trans tracks the current position
        [new,param] = translate (rotated, param.trans+d, param);
        param.trans = param.trans + d;
        str1 = print_param(param);
        % measure mutual_information
        if strcmp(param.myfunc_MI,'multiply')
            MI = mutual_information (pos1, index1, side1, new, side2, param, 0);
            MIt = mutual_information (pos1, index1, side1, new, side2, param, param.threshold_plot);
        elseif strcmp(param.myfunc_MI,'multiply_sqrt')
            MI = mutual_information_sqrt (pos1, index1, side1, new, side2, param, 0);
            MIt = mutual_information_sqrt (pos1, index1, side1, new, side2, param, param.threshold_plot);
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
            param.MItvec = [param.MItvec MIt];
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
                param.MItvec = [param.MItvec MIt];
                last_MI = MI;
            else
                % reject move
                str3 = sprintf('Reject %.3g > %.3g',rnd,p);
                param.rot = param.rot - r;
                param.trans = param.trans - d;
                tmp = param.MIvec(end);
                param.MIvec = [param.MIvec tmp];
                tmpt = param.MItvec(end);
                param.MItvec = [param.MItvec tmpt];
                
            end
        end
        pos_changes = pos_changes-1;
        str4 = sprintf('%d %d',Tchanges,pos_changes);
        str6 = sprintf('Final MI = %.3g',last_MI);
        str5 = print_param(param);
        val7 = param.trans - param.centroid;
        str7 = sprintf('final offset = [%6.6g0 %6.6g0 %6.6g0]',val7(1),val7(2),val7(3));
        w = find(param.centers>last_MI,1); % keep only first instance
        if ~isempty(w)
            str8 = sprintf('MI frac = %.3g,',param.cdf(w));
            param.cdfvec = [param.cdfvec param.cdf(w)];
        else
            str8 = 'MI frac = 1.0';
            param.cdfvec = [param.cdfvec 1];
        end
        dif = last_MI - max(param.nullMIvec);
        % last_MI > max(param.nullMIvec) =>dif>0
        % last_MI < max(param.nullMIvec) =>dif<0
        if dif > 0
            str9 = sprintf('siman exceeds null by %.3f',dif);
        else
            str9 = sprintf('null exceeds siman by %.3f',-dif);
        end
        disp(sprintf('%7s%75s  %84s  %40s  %22s  %84s  %20s %40s %20s %20s',str4,str0,str1,str2,str3,str5,str6,str7,str8,str9));
        param.Pvec = [param.Pvec p];
        param.transvec = [param.transvec; param.trans];
        param.rotvec = [param.rotvec; param.rot];
%         profile off
%         profile viewer
%         keyboard
    end
    %profile viewer;
    %keyboard
    T = lowerT(T,param);
    Tchanges = Tchanges-1;
    p = p * param.prate;
end
[rotated,param] = rotate (canonical,param.rot, param);
[new,param] = translate (rotated, param.trans, param);
end

function save_plot (h, param,f)
prefix = sprintf('%s%s_',param.savePath,param.timestamp);
str=sprintf('%s_%s.png',prefix,f);
fr=getframe(h);
[X, map] = frame2im(fr);
imwrite(X, str);
end

function save_stats (param)
prefix = sprintf('%s%s_',param.savePath,param.timestamp);
% plot Pvec
h = figure;
plot([1:numel(param.Pvec)],param.Pvec);
xlabel('iteration');
ylabel('probability of allowing a decrease in MI');
title('simulated annealing schedule');
if 0>1
    fname = sprintf('%s_p.fig',prefix)
    savefig(h,fname);
else
    str=sprintf('%s_p.png',prefix);
    f=getframe(gcf);
    [X, map] = frame2im(f);
    imwrite(X, str);
end

h = figure;
plot([1:numel(param.MIvec)],log10(param.MIvec));
xlabel('iteration');
ylabel('log (mutual information)');
title('evolution of mutual information during simulated annealing');
if 0>1
    fname = sprintf('%s_MI.fig',prefix)
    savefig(h,fname);
else
    str=sprintf('%s_MI.png',prefix);
    f=getframe(gcf);
    [X, map] = frame2im(f);
    imwrite(X, str);
end

f = figure;
plot([1:numel(param.cdfvec)],param.cdfvec);
xlabel('iteration');
ylabel('fraction of null distribution');
title('fraction of null distribution less than current mutual informaton');
if 0>1
    fname = sprintf('%s_frac_null.fig',prefix)
    savefig(h,fname);
else
    str=sprintf('%s_frac_null.png',prefix);
    ff=getframe(f);
    [X, map] = frame2im(ff);
    imwrite(X, str);
end

h = figure;
a = size(param.transvec);
subplot(1,3,1);
plot([1:a(1)],param.transvec(:,1));
xlabel('iteration');
ylabel('translation in dim 1 (um)');
subplot(1,3,2);
plot([1:a(1)],param.transvec(:,2));
xlabel('iteration');
ylabel('translation in dim 2 (um)');
title('trajectory of simulated annealing');
subplot(1,3,3);
plot([1:a(1)],param.transvec(:,3));
xlabel('iteration');
ylabel('translation in dim 3 (um)');
if 0>1
    fname = sprintf('%s_trans.fig',prefix)
    savefig(h,fname);
else
    %fname = sprintf('%s%s_p.fig',param.savePath,param.timestamp)
    str=sprintf('%s_trans.png',prefix);
    f=getframe(gcf);
    [X, map] = frame2im(f);
    imwrite(X, str);
end

h = figure;
a = size(param.rotvec);
subplot(1,3,1);
plot([1:a(1)],param.rotvec(:,1));
xlabel('iteration');
ylabel('rotation around dim 1 (radians)');
subplot(1,3,2);
plot([1:a(1)],param.rotvec(:,2));
xlabel('iteration');
ylabel('rotation around dim 2 (radians)');
title('trajectory of simulated annealing');
subplot(1,3,3);
plot([1:a(1)],param.rotvec(:,3));
xlabel('iteration');
ylabel('rotation around dim 3 (radians)');
if 0>1
    fname = sprintf('%s_rot.fig',prefix)
    savefig(h,fname);
else
    %fname = sprintf('%s%s_p.fig',param.savePath,param.timestamp)
    str=sprintf('%s_rot.png',prefix);
    f=getframe(gcf);
    [X, map] = frame2im(f);
    imwrite(X, str);
end


end


function out = print_param (param)
a=param.trans;
b=param.rot;
out = sprintf('trans = [%6.6g0 %6.6g0 %6.6g0], rot = [%7.6g0 %7.6g0 %7.6g0]',a(1),a(2),a(3),b(1),b(2),b(3));
end


function mi = mutual_information (pos1, index1, side1, new, side2, param, threshold)
mi = 0;
s = size(side1);
scale = [1/param.voxel_y 0 0 ; 0 1/param.voxel_x 0 ; 0 0 1/param.voxel_z];
tmp = ceil(new*scale);
for i=1:length(tmp)
    % convert position to a,b,c
    a = tmp(i,1);
    b = tmp(i,2);
    c = tmp(i,3);
    % check
    if a>s(1) || b>s(2) || c>s(3) || a<1 || b<1 || c<1
        continue;
    else
        i1 = side1(a,b,c);
        i2 = side2(param.index2(i));
        if threshold
            if i1<threshold && i2<threshold
                continue;
            end
        end
        mi = mi + double(double(i1)*double(i2));
    end
end
end

function mi = mutual_information_sqrt (pos1, index1, side1, new, side2, param, threshold)
mi = 0;
s = size(side1);
scale = [1/param.voxel_y 0 0 ; 0 1/param.voxel_x 0 ; 0 0 1/param.voxel_z];
tmp = ceil(new*scale);
for i=1:length(tmp)
    % convert position to a,b,c
    a = tmp(i,1);
    b = tmp(i,2);
    c = tmp(i,3);
    % check
    if a>s(1) || b>s(2) || c>s(3) || a<1 || b<1 || c<1
        continue;
    else
        i1 = side1(a,b,c);
        i2 = side2(param.index2(i));
        if threshold
            if i1<threshold && i2<threshold
                continue;
            end
        end
        mi = mi + sqrt(double(double(i1)*double(i2)));
    end
end
end

function plot_pos (pos1, pos2, side1, side2, new, param, f)

glyph = {'.','.'};
colors = {[0 0 1],
    [0.8 0.8 0.8],
    [1 0.8 0.8],};

N = 1000;

%index1 = find(side1>param.threshold);
%index2 = find(side2>param.threshold);

h = figure;
subplot(2,2,2);
hold on;

n = ceil(length(param.index2)/N);
for i=1:n:length(param.index2)
    val = 1-single(side2(param.index2(i)))/2^16;
    plot(pos2(i,2), -pos2(i,1), glyph{2},'Color',[val val val]);
end

n = ceil(length(param.index2)/N);
for i=1:n:length(param.index2)
    val = 1-single(side2(param.index2(i)))/2^16;
    plot(new(i,2), -new(i,1), glyph{2},'Color',[1 val val]);
end

n = ceil(length(param.index1)/N);
for i=1:n:length(param.index1)
    val = 1-single(side1(param.index1(i)))/2^16;
    plot(pos1(i,2), -pos1(i,1), glyph{2},'Color',[val val 1]);
end

% text(300,-250,'LFM1','FontSize',12,'Color',[0 0 1]);
% text(300,-270,'LFM2','FontSize',12,'Color',[0.8 0.8 0.8]);
% text(300,-290,'DLFM','FontSize',12,'Color',[1 0 0]);

q = param.centroid; % centroid2
plot( q(2), -q(1), 'o','Color',colors{2});
hold off;
xlabel('two [um]');
ylabel('one [um]');
axis equal;


subplot(2,2,1);
hold on;

n = ceil(length(param.index2)/N);
for i=1:n:length(param.index2)
    val = 1-single(side2(param.index2(i)))/2^16;
    plot(pos2(i,3), -pos2(i,1), glyph{2},'Color',[val val val]);
end

n = ceil(length(param.index2)/N);
for i=1:n:length(param.index2)
    val = 1-single(side2(param.index2(i)))/2^16;
    plot(new(i,3), -new(i,1), glyph{2},'Color',[1 val val]);
end

n = ceil(length(param.index1)/N);
for i=1:n:length(param.index1)
    val = 1-single(side1(param.index1(i)))/2^16;
    plot(pos1(i,3), -pos1(i,1), glyph{2},'Color',[val val 1]);
end

% text(300,-250,'LFM1','FontSize',12,'Color',[0 0 1]);
% text(300,-270,'LFM2','FontSize',12,'Color',[0.8 0.8 0.8]);
% text(300,-290,'DLFM','FontSize',12,'Color',[1 0 0]);

q = param.centroid; % centroid2
plot( q(3), -q(1), 'o','Color',colors{2});
hold off;
axis equal;
xlabel('three [um]');
ylabel('one [um]');

subplot(2,2,4);
hold on;

n = ceil(length(param.index2)/N);
for i=1:n:length(param.index2)
    val = 1-single(side2(param.index2(i)))/2^16;
    plot(pos2(i,2), pos2(i,3), glyph{2},'Color',[val val val]);
end

n = ceil(length(param.index2)/N);
for i=1:n:length(param.index2)
    val = 1-single(side2(param.index2(i)))/2^16;
    plot(new(i,2), new(i,3), glyph{2},'Color',[1 val val]);
end

n = ceil(length(param.index1)/N);
for i=1:n:length(param.index1)
    val = 1-single(side1(param.index1(i)))/2^16;
    plot(pos1(i,2), pos1(i,3), glyph{2},'Color',[val val 1]);
end


text(-800,400,'LFM1','FontSize',12,'Color',[0 0 1]);
text(-800,300,'LFM2','FontSize',12,'Color',[0.8 0.8 0.8]);
text(-800,200,'DLFM','FontSize',12,'Color',[1 0 0]);

q = param.centroid; % centroid2
plot( q(2), q(3), 'o','Color',colors{2});
hold off;

axis equal;
xlabel('two [um]');
ylabel('three [um]');

prefix = sprintf('%s%s_',param.savePath,param.timestamp);
str=sprintf('%s_%s.png',prefix,f);
fr=getframe(h);
[X, map] = frame2im(fr);
imwrite(X, str);
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


function i = intensity (index, side)
[a,b,c] = ind2sub(size(side),index);
i = double(side(a,b,c));
if i==0
    disp('WTF!');
    keyboard;
end
end


function [out,param] = rotate (pos, angle, param)
rot = single( rotation_matrix (angle) );
out = pos*rot;
end


function [out,param] = translate (pos, delta, param)
D = ones(size(pos),'single').*delta;
out = pos + D;
end

function p = calc_boltzman_p (MI, T, param)
val = -param.scale/T;
if val > 0.0
    p = 1.0;
else
    p = exp(val);
end
if p<0.0 || p>1.0
    disp('WTF?');
    keyboard
end
end


function out = calc_centroid (side, param)
s = size(side);
out = s/2.*[param.voxel_y param.voxel_x param.voxel_z];
end


function pos = init_pos (linear_i,side, param)
[a,b,c] = ind2sub(size(side),linear_i);
y = single( (a-0.5) * param.voxel_y );
x = single( (b-0.5) * param.voxel_x );
z = single( (c-0.5) * param.voxel_z );
pos = [y x z];
end


function L = linear_new (new, side1, param)
L = zeros([length(new) 1]);
for i=1:length(new)
    % convert position to a,b,c
    y = new(i,1);
    a = ceil(y/param.voxel_y);
    x = new(i,2);
    b = ceil(x/param.voxel_x);
    z = new(i,3);
    c = ceil(z/param.voxel_z);
    % check
    s = size(side1);
    if a>s(1) || b>s(2) || c>s(3) || a<1 || b<1 || c<1
        L(i) = -1;
    else
        % convert a,b,c to linear index
        ind = sub2ind(size(side1),a,b,c);
        L(i) = ind;
    end
end
end


function save_vol (A, savePath, outFile)
imwrite( squeeze(A(:,:,1)), [savePath outFile]);
for k = 2:size(A,3),
    imwrite(squeeze(A(:,:,k)),  [savePath outFile], 'WriteMode', 'append');
end
end


function [d,r] = perturb (param, gain)
d1 = random('unif',-1,1) * param.trans_amp * gain;
d2 = random('unif',-1,1) * param.trans_amp * gain;
d3 = random('unif',-1,1) * param.trans_amp * gain;
r1 = random('unif',-1,1) * param.rot_amp * gain;
r2 = random('unif',-1,1) * param.rot_amp * gain;
r3 = random('unif',-1,1) * param.rot_amp * gain;
d = [d1 d2 d3];
%r = [r1 r2 r3];
r = [r1 0 0];
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

function param = reeset (param)
param.voxel_x = 0.323; % um
param.voxel_y = 0.323; % um
param.voxel_z = 4.0; % um
param.interp = 8;
param.myfunc_combine = 'multiply_sqrt';
param.myfunc_MI = 'multiply_sqrt';

param.prate = -1;
param.Trate = 1e-1;
param.final_p = 1e-2;
param.init_p = 0.3;
param.scale_trans = 1;
param.scale_rot = 1;
param.trans_amp = param.scale_trans * param.voxel_x; % um
param.rot_amp = param.scale_rot * pi/800; % radians
%param.trans_amp = 1.0; % um
%param.rot_amp = pi/80; % radians

% amount to clip from periphery in pixels
% dim1 min, dim1 max, dim2 min, dim2 max, dim3 min, dim3 max,
param.clip = [0,0,0,0,0,0];

param.power = 1;
param.T0 = -1;
param.TC0 = 30;
param.MC0 = 30;
%param.TC0 = 1;
%param.MC0 = 900;

param.psf   = [-1.0 -1.0 -1.0]; %um
%param.offset = [-8 -30 0];
param.offset = [-5 0 -60];
param.trans = [0 0 0]; % um
%param.angle   = [-1.2*pi/2 0 0]; % radians
param.angle   = [-1.0*pi/2 0 0]; % radians
param.rot   = [0 0 0]; % radians

param.lfdisplay = false;

%param.centroid = [169.8834 184.3252 129.1565];
param.centroid = [-1 -1 -1];

% indices of side2 that are being tracked
%param.threshold = 20;
param.threshold1 = 40;
param.threshold2 = 40;
param.threshold_plot = 40;
param.index2 = [];

param.rapid = false;
param.savevol = true;
param.plot = true;
%param.savePath      = '/Users/justin/Desktop/LFM volume registration/';
param.savePath      = '/home/jkinney/Desktop/LFM volume registration/registration/param_sweep/';
%param.savePath      = '/tmp/';
param.inputFilePath = '/home/jkinney/Desktop/LFM volume registration/from_nikita/';
param.inputFileName = {
    'raw_side1.mat',
    'raw_side2.mat'};
end


function print_fraction (index, side, str)
i = 0;
for j=1:length(index)
    i = i + single(side(index(j)));
end
iT = single(sum(sum(sum(side))));

fi = i/iT;
fv = single(length(index))/single(numel(side));
disp(sprintf('%s: %d of %d voxels (%7.3g), %d of %d total intensity (%7.3g)',...
    str,length(index),numel(side),fv,...
    i,iT,fi));
end

function save_1d_summed_projections (side1, side2, pos1, pos2, new, param)
glyph = {'.','.'};
colors = {[1 0 0],
    [0.8 0.8 0.8],
    [0 0 1],
    [0.4 0.4 0.4],};

onetwo = squeeze(sum(side1,3));
twothree = squeeze(sum(side1,1));
side1_d1 = squeeze(sum(onetwo,2));
side1_d2 = squeeze(sum(onetwo,1));
side1_d3 = squeeze(sum(twothree,1));

onetwo = squeeze(sum(side2,3));
twothree = squeeze(sum(side2,1));
side2_d1 = squeeze(sum(onetwo,2));
side2_d2 = squeeze(sum(onetwo,1));
side2_d3 = squeeze(sum(twothree,1));


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
        onetwo(abc_new(i,1),abc_new(i,2)) = onetwo(abc_new(i,1),abc_new(i,2)) + double(side2(param.index2(i)));
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(new)
    if abc_new(i,2) > 0 && abc_new(i,3) > 0
        twothree(abc_new(i,2),abc_new(i,3)) = twothree(abc_new(i,2),abc_new(i,3)) + double(side2(param.index2(i)));
    end
end
% convert 2d projections to 1d
new_d1 = squeeze(sum(onetwo,2));
new_d2 = squeeze(sum(onetwo,1));
new_d3 = squeeze(sum(twothree,1));

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
        onetwo(abc_pos1(i,1),abc_pos1(i,2)) = onetwo(abc_pos1(i,1),abc_pos1(i,2)) + double(side1(param.index1(i)));
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(pos1)
    if abc_pos1(i,2) > 0 && abc_pos1(i,3) > 0
        twothree(abc_pos1(i,2),abc_pos1(i,3)) = twothree(abc_pos1(i,2),abc_pos1(i,3)) + double(side1(param.index1(i)));
    end
end
pos1_d1 = squeeze(sum(onetwo,2));
pos1_d2 = squeeze(sum(onetwo,1));
pos1_d3 = squeeze(sum(twothree,1));

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
        onetwo(abc_pos2(i,1),abc_pos2(i,2)) = onetwo(abc_pos2(i,1),abc_pos2(i,2)) + double(side2(param.index2(i)));
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(pos2)
    if abc_pos2(i,2) > 0 && abc_pos2(i,3) > 0
        twothree(abc_pos2(i,2),abc_pos2(i,3)) = twothree(abc_pos2(i,2),abc_pos2(i,3)) + double(side2(param.index2(i)));
    end
end
pos2_d1 = squeeze(sum(onetwo,2));
pos2_d2 = squeeze(sum(onetwo,1));
pos2_d3 = squeeze(sum(twothree,1));

if 0 < 1
    T = 'Summed intensity projection (scaled by summed intensity)';
    scale_d1 = [sum(pos1_d1) sum(pos2_d1) sum(side1_d1) sum(side2_d1) sum(new_d1)];
    scale_d2 = [sum(pos1_d2) sum(pos2_d2) sum(side1_d2) sum(side2_d2) sum(new_d2)];
    scale_d3 = [sum(pos1_d3) sum(pos2_d3) sum(side1_d3) sum(side2_d3) sum(new_d3)];
    %disp(T);
elseif 0 > 1
    T = 'Summed intensity projection (scaled by max intensity)';
    scale_d1 = [max(pos1_d1) max(pos2_d1) max(side1_d1) max(side2_d1) max(new_d1)];
    scale_d2 = [max(pos1_d2) max(pos2_d2) max(side1_d2) max(side2_d2) max(new_d2)];
    scale_d3 = [max(pos1_d3) max(pos2_d3) max(side1_d3) max(side2_d3) max(new_d3)];
    %disp(T);
else
    T = 'Summed intensity projection (not scaled)';
    scale_d1 = [1 1 1 1 1];
    scale_d2 = [1 1 1 1 1];
    scale_d3 = [1 1 1 1 1];
    %disp(T);
end



f = figure;
set(gcf,'Position',[79          18        1270         940]);
subplot(3,1,1);
plot(pos1_d1/scale_d1(1),'o','Color',colors{1});
hold on;
plot(pos2_d1/scale_d1(2),'o','Color',colors{4});
plot(side1_d1/scale_d1(3),'.','Color',colors{2});
plot(side2_d1/scale_d1(4),'*','Color',colors{2});
plot(new_d1/scale_d1(5),'-','Color',colors{4});
xlabel('first dimension [pixels]');
ylabel('intensity');
hold off;
legend('pos1','pos2','side1','side2','new');
title(T);
%xlim([100 1100]);

subplot(3,1,2);
plot(pos1_d2/scale_d2(1),'o','Color',colors{1});
hold on;
plot(pos2_d2/scale_d2(2),'o','Color',colors{4});
plot(side1_d2/scale_d2(3),'.','Color',colors{2});
plot(side2_d2/scale_d2(4),'*','Color',colors{2});
plot(new_d2/scale_d2(5),'-','Color',colors{4});
xlabel('second dimension [pixels]');
ylabel('intensity');
hold off;
legend('pos1','pos2','side1','side2','new');
%xlim([100 1100]);
%ylim([0 0.005]);

subplot(3,1,3);
plot(pos1_d3/scale_d3(1),'o','Color',colors{1});
hold on;
plot(pos2_d3/scale_d3(2),'o','Color',colors{4});
plot(side1_d3/scale_d3(3),'.','Color',colors{2});
plot(side2_d3/scale_d3(4),'*','Color',colors{2});
plot(new_d3/scale_d3(5),'-','Color',colors{4});
xlabel('third dimension [pixels]');
ylabel('intensity');
hold off;
legend('pos1','pos2','side1','side2','new');
%xlim([100 1100]);
%ylim([0 0.005]);

str=sprintf('%s%s_summed_intensity_projection.png',param.savePath,param.timestamp);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);
end

function save_1d_max_projections (side1, side2, pos1, pos2, new, param, str)
glyph = {'.','.'};
colors = {[0 0 1],
    [0.8 0.8 0.8],
    [0 0 1],
    [1 0 0],};

onetwo = squeeze(max(side1,[],3));
twothree = squeeze(max(side1,[],1));
side1_d1 = single(squeeze(max(onetwo,[],2)));
side1_d2 = single(squeeze(max(onetwo,[],1)));
side1_d3 = single(squeeze(max(twothree,[],1)));

onetwo = squeeze(max(side2,[],3));
twothree = squeeze(max(side2,[],1));
side2_d1 = single(squeeze(max(onetwo,[],2)));
side2_d2 = single(squeeze(max(onetwo,[],1)));
side2_d3 = single(squeeze(max(twothree,[],1)));

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
        onetwo(abc_new(i,1),abc_new(i,2)) = max([onetwo(abc_new(i,1),abc_new(i,2)) double(side2(param.index2(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(new)
    if abc_new(i,2) > 0 && abc_new(i,3) > 0
        twothree(abc_new(i,2),abc_new(i,3)) = max([twothree(abc_new(i,2),abc_new(i,3))  double(side2(param.index2(i)))]);
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
        onetwo(abc_pos1(i,1),abc_pos1(i,2)) = max([onetwo(abc_pos1(i,1),abc_pos1(i,2)) double(side1(param.index1(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(pos1)
    if abc_pos1(i,2) > 0 && abc_pos1(i,3) > 0
        twothree(abc_pos1(i,2),abc_pos1(i,3)) = max([twothree(abc_pos1(i,2),abc_pos1(i,3)) double(side1(param.index1(i)))]);
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
        onetwo(abc_pos2(i,1),abc_pos2(i,2)) = max([onetwo(abc_pos2(i,1),abc_pos2(i,2)) double(side2(param.index2(i)))]);
    end
end
twothree = zeros(ns(2),ns(3));
for i=1:length(pos2)
    if abc_pos2(i,2) > 0 && abc_pos2(i,3) > 0
        twothree(abc_pos2(i,2),abc_pos2(i,3)) = max([twothree(abc_pos2(i,2),abc_pos2(i,3)) double(side2(param.index2(i)))]);
    end
end
pos2_d1 = squeeze(max(onetwo,[],2));
pos2_d2 = squeeze(max(onetwo,[],1));
pos2_d3 = squeeze(max(twothree,[],1));

if 0 < 1
    T = '      (scaled by total intensity in sample)';
    scale_d1 = [sum(pos1_d1) sum(pos2_d1) sum(side1_d1) sum(side2_d1) sum(new_d1)];
    scale_d2 = [sum(pos1_d2) sum(pos2_d2) sum(side1_d2) sum(side2_d2) sum(new_d2)];
    scale_d3 = [sum(pos1_d3) sum(pos2_d3) sum(side1_d3) sum(side2_d3) sum(new_d3)];
    yl = 'normalized intensity';
    %disp(T);
elseif 0 < 1
    T = '      (scaled by max intensity in each projected volume)';
    scale_d1 = single([max(pos1_d1) max(pos2_d1) max(side1_d1) max(side2_d1) max(new_d1)]);
    scale_d2 = single([max(pos1_d2) max(pos2_d2) max(side1_d2) max(side2_d2) max(new_d2)]);
    scale_d3 = single([max(pos1_d3) max(pos2_d3) max(side1_d3) max(side2_d3) max(new_d3)]);
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
semilogy(side1_d1/scale_d1(3),'.','Color',colors{1});
hold on;
semilogy(side2_d1/scale_d1(4),'.','Color',colors{2});
%semilogy(pos1_d1/scale_d1(1),'o','Color',colors{1});
%semilogy(pos2_d1/scale_d1(2),'o','Color',colors{2});
semilogy(new_d1/scale_d1(5),'-','Color',colors{4});
xlabel('first dimension [pixels]');
ylabel(yl);
hold off;
%legend('side1','side2','pos1','pos2','new');
%legend('pos1','pos2','new');
legend('LFM1','LFM2','LFM2 coarse reg');
title([str T],'Interpreter','none');
%xlim([100 1100]);

subplot(3,1,2);
semilogy(side1_d2/scale_d2(3),'.','Color',colors{1});
hold on;
semilogy(side2_d2/scale_d2(4),'.','Color',colors{2});
%semilogy(pos1_d2/scale_d2(1),'o','Color',colors{1});
%semilogy(pos2_d2/scale_d2(2),'o','Color',colors{2});
semilogy(new_d2/scale_d2(5),'-','Color',colors{4});
xlabel('second dimension [pixels]');
ylabel(yl);
hold off;
%legend('side1','side2','pos1','pos2','new');
%legend('pos1','pos2','new');
legend('LFM1','LFM2','LFM2 coarse reg');
%xlim([100 1100]);
%ylim([0 0.005]);

subplot(3,1,3);
semilogy(side1_d3/scale_d3(3),'.','Color',colors{1});
hold on;
semilogy(side2_d3/scale_d3(4),'.','Color',colors{2});
%semilogy(pos1_d3/scale_d3(1),'o','Color',colors{1});
%semilogy(pos2_d3/scale_d3(2),'o','Color',colors{2});
semilogy(new_d3/scale_d3(5),'-','Color',colors{4});
xlabel('third dimension [pixels]');
ylabel(yl);
hold off;
%legend('side1','side2','pos1','pos2','new');
%legend('pos1','pos2','new');
legend('LFM1','LFM2','LFM2 coarse reg');
%xlim([100 1100]);
%ylim([0 0.005]);

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);
end


function save_2d_contour_plots (side1, side2, pos1, pos2, new, param, str)
glyph = {'.','.'};
colors = {[1 0 0],
    [0.8 0.8 0.8],
    [0 0 1],
    [0.4 0.4 0.4],};

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
            double(side2(param.index2(i)))]);
    end
end
onethree_new = zeros(ns(1),ns(3));
for i=1:length(new)
    if abc_new(i,1) > 0 && abc_new(i,3) > 0
        onethree_new(abc_new(i,1),abc_new(i,3)) ...
            = max([onethree_new(abc_new(i,1),abc_new(i,3))...
            double(side2(param.index2(i)))]);
    end
end
threetwo_new = zeros(ns(3),ns(2));
for i=1:length(new)
    if abc_new(i,3) > 0 && abc_new(i,2) > 0
        threetwo_new(abc_new(i,3),abc_new(i,2))...
            = max([threetwo_new(abc_new(i,3),abc_new(i,2))...
            double(side2(param.index2(i)))]);
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
            double(side1(param.index1(i)))]);
    end
end
onethree_pos1 = zeros(ns(1),ns(3));
for i=1:length(pos1)
    if abc_pos1(i,1) > 0 && abc_pos1(i,3) > 0
        onethree_pos1(abc_pos1(i,1),abc_pos1(i,3))...
            = max([onethree_pos1(abc_pos1(i,1),abc_pos1(i,3))...
            double(side1(param.index1(i)))]);
    end
end
threetwo_pos1 = zeros(ns(3),ns(2));
for i=1:length(pos1)
    if abc_pos1(i,3) > 0 && abc_pos1(i,2) > 0
        threetwo_pos1(abc_pos1(i,3),abc_pos1(i,2))...
            = max([threetwo_pos1(abc_pos1(i,3),abc_pos1(i,2))...
            double(side1(param.index1(i)))]);
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
            double(side2(param.index2(i)))]);
    end
end
onethree_pos2 = zeros(ns(1),ns(3));
for i=1:length(pos2)
    if abc_pos2(i,2) > 0 && abc_pos2(i,3) > 0
        onethree_pos2(abc_pos2(i,1),abc_pos2(i,3)) ...
            = max([onethree_pos2(abc_pos2(i,1),abc_pos2(i,3)) ...
            double(side2(param.index2(i)))]);
    end
end
threetwo_pos2 = zeros(ns(3),ns(2));
for i=1:length(pos2)
    if abc_pos2(i,3) > 0 && abc_pos2(i,2) > 0
        threetwo_pos2(abc_pos2(i,3),abc_pos2(i,2)) ...
            = max([threetwo_pos2(abc_pos2(i,3),abc_pos2(i,2)) ...
            double(side2(param.index2(i)))]);
    end
end

f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
hold on;
a = size(onethree_new);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onethree_new,v,'Color',colors{1});
a = size(onethree_pos1);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,onethree_pos1,v,'Color',colors{3});
a = size(onethree_pos2);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onethree_pos2,v,'Color',colors{2});
xlabel('dim three [pixels]');
ylabel('dim one [pixels]');
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

subplot(2,2,2);
hold on;
a = size(onetwo_new);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onetwo_new,v,'Color',colors{1});
a = size(onetwo_pos1);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,onetwo_pos1,v,'Color',colors{3});
a = size(onetwo_pos2);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,onetwo_pos2,v,'Color',colors{2});
xlabel('dim two [pixels]');
ylabel('dim one [pixels]');
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

subplot(2,2,4);
hold on;
a = size(threetwo_new);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,threetwo_new,v,'Color',colors{1});
a = size(threetwo_pos1);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int1;
v = [ clevel, clevel ];
contour(X,Y,threetwo_pos1,v,'Color',colors{3});
a = size(threetwo_pos2);
x = [1:a(2)];
y = [1:a(1)];
[X,Y] = meshgrid(x,y);
clevel = param.contour_int2;
v = [ clevel, clevel ];
contour(X,Y,threetwo_pos2,v,'Color',colors{2});
xlabel('dim two [pixels]');
ylabel('dim three [pixels]');
daspect([1,1,1]);
hold off;
set(gca,'Ydir','reverse');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.25,0.4,'LFM1','FontSize',12,'Color',colors{3},'Interpreter','none');
text(0.25,0.35,'LFM2','FontSize',12,'Color',colors{2},'Interpreter','none');
text(0.25,0.3,'LFM2 coarse reg','FontSize',12,'Color',colors{1},'Interpreter','none');
text(0.5,0.92,str,'FontSize',12,'Color',[0 0 0] ,'Interpreter','none');

text(0.5,0.8,sprintf('contour at intensity level %4.0f',param.contour_int1),'FontSize',12,'Color',colors{3},'Interpreter','none');
text(0.5,0.75,sprintf('contour at intensity level %4.0f',param.contour_int2),'FontSize',12,'Color',colors{2},'Interpreter','none');
text(0.5,0.7,sprintf('contour at intensity level %4.0f',param.contour_int2),'FontSize',12,'Color',colors{1},'Interpreter','none');

str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);
end


% function save_2d_projections (side1, side2, pos1, pos2, new, param, flag, str)
% %f = figure();
% %set(gcf,'Position',[1367          19        1184         952]);
% f = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% yz = squeeze(sum(side1,2));
% subplot(2,6,1);
% dr = ceil(log2(max(max(yz))));
% imagesc(yz,[0 2^dr]);
% xlabel('three [pixels]');
% ylabel('one [pixels]');
% title('LFM1');
% colorbar();
% daspect([1,1,1]);
% 
% xy = squeeze(sum(side1,3));
% subplot(2,6,2);
% dr = ceil(log2(max(max(xy))));
% imagesc(xy,[0 2^dr]);
% xlabel('two [pixels]');
% ylabel('one [pixels]');
% title('LFM1');
% colorbar();
% daspect([1,1,1]);
% 
% xz = squeeze(sum(side1,1))';
% subplot(2,6,8);
% dr = ceil(log2(max(max(xz))));
% imagesc(xz,[0 2^dr]);
% xlabel('two [pixels]');
% ylabel('three [pixels]');
% title('LFM1');
% colorbar();
% daspect([1,1,1]);
% 
% 
% yz = squeeze(sum(side2,2));
% subplot(2,6,3);
% dr = ceil(log2(max(max(yz))));
% imagesc(yz,[0 2^dr]);
% xlabel('three [pixels]');
% ylabel('one [pixels]');
% title('LFM2');
% colorbar();
% daspect([1,1,1]);
% 
% xy = squeeze(sum(side2,3));
% subplot(2,6,4);
% dr = ceil(log2(max(max(xy))));
% imagesc(xy,[0 2^dr]);
% xlabel('two [pixels]');
% ylabel('one [pixels]');
% title('LFM2');
% colorbar();
% daspect([1,1,1]);
% 
% xz = squeeze(sum(side2,1))';
% subplot(2,6,10);
% dr = ceil(log2(max(max(xz))));
% imagesc(xz,[0 2^dr]);
% xlabel('two [pixels]');
% ylabel('three [pixels]');
% title('LFM2');
% colorbar();
% daspect([1,1,1]);
% 
% if flag
%     yz = squeeze(sum(new,2));
%     subplot(2,6,5);
%     dr = ceil(log2(max(max(yz))));
%     imagesc(yz,[0 2^dr]);
%     xlabel('three [pixels]');
%     ylabel('one [pixels]');
%     title('DLFM');
%     colorbar();
%     daspect([1,1,1]);
%     
%     xy = squeeze(sum(new,3));
%     subplot(2,6,6);
%     dr = ceil(log2(max(max(xy))));
%     imagesc(xy,[0 2^dr]);
%     xlabel('two [pixels]');
%     ylabel('one [pixels]');
%     title('DLFM');
%     colorbar();
%     daspect([1,1,1]);
%     
%     xz = squeeze(sum(new,1))';
%     subplot(2,6,12);
%     dr = ceil(log2(max(max(xz))));
%     imagesc(xz,[0 2^dr]);
%     xlabel('two [pixels]');
%     ylabel('three [pixels]');
%     title('DLFM');
%     colorbar();
%     daspect([1,1,1]);
% end
% 
% 
% % set all colorbars to same range
% crange = [];
% i = [1 2 3 4 8 10];
% if flag
%     i = [i 5 6 12];
% end
% 
% for j=i
%     subplot(2,6,j);
%     crange = [crange;caxis];
% end
% m = max(crange);
% for j=i
%     subplot(2,6,j);
%     caxis(m);
% end
% 
% % save figure
% str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
% ff=getframe(f);
% [X, map] = frame2im(ff);
% imwrite(X, str);
% end

function save_2d_max_projections (side1, side2, new, param, flag, str)

f = figure('units','normalized','outerposition',[0 0 1 1]);
yz = squeeze(max(side1,[],2));
subplot(2,6,1);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM1');
%colorbar();
daspect([1,1,1]);

xy = squeeze(max(side1,[],3));
subplot(2,6,2);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM1');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(side1,[],1))';
subplot(2,6,8);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
title('LFM1');
colorbar();
daspect([1,1,1]);

yz = squeeze(max(side2,[],2));
subplot(2,6,3);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
title('LFM2');
%colorbar();
daspect([1,1,1]);

xy = squeeze(max(side2,[],3));
subplot(2,6,4);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
title('LFM2');
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(side2,[],1))';
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

% % set all colorbars to same range
% crange = [];
% i = [1 2 3 4 8 10];
% if flag
%     i = [i 5 6 12];
% end
% 
% for j=i
%     subplot(2,6,j);
%     crange = [crange;caxis];
% end
% m = max(crange);
% for j=i
%     subplot(2,6,j);
%     caxis(m);
% end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.15,0.97,[param.inputFilePath1 param.inputFileName],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
text(0.4,0.97,[param.inputFilePath2 param.inputFileName],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
if flag
    text(0.65,0.97,[param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
end
text(0.4,0.2,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);
end


function save_2d_max_projections_compact (side1, side2, new, param, flag, str)

f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1);
yz = squeeze(max(side1,[],2));
xy = squeeze(max(side1,[],3));
xz = squeeze(max(side1,[],1))';
a = size(side1);
big_image = [a(3)+a(1) a(3)+a(2)];
imagesc(zeros(big_image));
a = size(yz);
x = [1 a(2)];
y = [1 a(1)];
imagesc('XData',x,'YData',y,'CData',yz);
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
yz = squeeze(max(side2,[],2));
xy = squeeze(max(side2,[],3));
xz = squeeze(max(side2,[],1))';
a = size(side2);
big_image = [a(3)+a(1) a(3)+a(2)];
imagesc(zeros(big_image));
a = size(yz);
x = [1 a(2)];
y = [1 a(1)];
imagesc('XData',x,'YData',y,'CData',yz);
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
    imagesc('XData',x,'YData',y,'CData',yz);
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
text(0.15,0.97,[param.inputFilePath1 param.inputFileName],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
text(0.4,0.97,[param.inputFilePath2 param.inputFileName],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
if flag
    text(0.65,0.97,[param.savePath param.timestamp],'FontSize',8,'Color',[0 0 0],'Interpreter','none');
end
text(0.4,0.2,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s%s_%s.png',param.savePath,param.timestamp, str);
ff=getframe(f);
[X, map] = frame2im(ff);
imwrite(X, str);
end



%% PARKING LOT
% For each temperature, run for long time to see equilibration
%
% for each threshold, report fraction of voxels and total intensity used
% to measure registration.
%
% Use uniform sampling of threshold to allow comparisons across thresholds.
