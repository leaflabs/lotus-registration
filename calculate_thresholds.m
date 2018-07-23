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
    
    %keyboard
    
    str=sprintf('%s%s_%s.png',param.savePath,param.timestamp,str);
    save_plot(f,str);
end
end