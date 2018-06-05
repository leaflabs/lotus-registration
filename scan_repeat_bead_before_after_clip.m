clear all;
close all;

%% load parameters
param = init_param();

[~,hostname] = system('hostname')

if ~isempty(strfind(hostname, 'Justins-Mac'))
	ppath = '/Users/justin/Desktop/DDLFM';
elseif ~isempty(strfind(hostname, 'willis'))
	ppath = '/home/jkinney/Desktop/DDLFM';
else
	disp('WTF?');
    keyboard
end

voxel_size = [0.323 0.323 4/8]; % um
ddlfm_path = 'bead/20180328/registration/';
opath = [ppath '/' ddlfm_path];
param.interp = 1;

savePath = [opath 'analysis/'];
param.savePath = savePath;
if ~exist(savePath,'dir')
    status = mkdir(savePath);
    if status == 1
        disp(['Created folder: ' savePath]);
    else
        disp(['Error attempting to create folder:' savePath]);
        status
        exit;
    end
end

timestamp = [datestr(datetime('now'),'yyyymmdd_HHMMSS')];
param.timestamp = timestamp;
fname = sprintf('%sscan_repeat_bead_before_after_clip_%s.log',savePath, timestamp);
diary(fname)
tic

% % row col z
% % y x z
% % dim 1 2 3
p = [opath '*sqrt.mat'];
d = dir(p);
fprintf('%d mat files found for\n%s\n\n',numel(d),p);
N = length(d);
for i=1:N
    fprintf('file %d of %d\n',i,N);
    f = [d(i).folder '/' d(i).name]
    % load volume
    fprintf('\nLoading volume %s\n\n',f);
    LFM1 = loadData(f, param);
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
    save_2d_max_projection(LFM1, param, '2d_max_projection',d(i).name);
end

diary off;

%% Functions

function save_2d_max_projection (LFM1, param, str, mytitle)

f = figure('units','normalized','outerposition',[0 0 1 1]);
yz = squeeze(max(LFM1,[],2));
subplot(2,2,1);
dr = ceil(log2(single(max(max(yz)))));
imagesc(yz,[0 2^dr]);
xlabel('three [pixels]');
ylabel('one [pixels]');
%colorbar();
daspect([1,1,1]);
set(gca,'XDir','reverse');

xy = squeeze(max(LFM1,[],3));
subplot(2,2,2);
dr = ceil(log2(single(max(max(xy)))));
imagesc(xy,[0 2^dr]);
xlabel('two [pixels]');
ylabel('one [pixels]');
%title(mytitle);
%colorbar();
daspect([1,1,1]);

xz = squeeze(max(LFM1,[],1))';
subplot(2,2,4);
dr = ceil(log2(single(max(max(xz)))));
imagesc(xz,[0 2^dr]);
xlabel('two [pixels]');
ylabel('three [pixels]');
colorbar();
daspect([1,1,1]);

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.1,0.98,mytitle,'FontSize',12,'Color',[0 0 0],'Interpreter','none');
%text(0.4,0.1,str,'FontSize',8,'Color',[0 0 0],'Interpreter','none');
drawnow

% save figure
str=sprintf('%s_%s_%s_%s.png',param.savePath,param.timestamp, mytitle(1:end-3),str);
print(f,str,'-dpng');
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

