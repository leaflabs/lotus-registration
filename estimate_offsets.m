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