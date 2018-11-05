function out = loadData (f, param, mystr)
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
    if 0>1
    out = interpolate (XguessSAVE1,param);
    else
        out = uint16(arrayInterp(single(XguessSAVE1), 3, 8, 'linear'));
    end
    if param.plot
        twothree = squeeze(sum(XguessSAVE1,1));
        XguessSAVE1_d3 = squeeze(sum(twothree,1));
        twothree = squeeze(sum(out,1));
        out_d3 = single(squeeze(sum(twothree,1)));
        f = figure;
        hold on;
        plot(out_d3,'o-');
        plot(linspace(1,numel(out_d3),numel(XguessSAVE1_d3)),XguessSAVE1_d3,'o-');
        hold off;
        legend('interpolated','original');
        xlabel('dim 3 [um]');
        ylabel('intensity [-]');
        title(['Summed intensity projection onto dim3 for ' mystr]);
        % add centroid values
        s = size(out);
        centroid_interp = s/2.*[param.voxel_y param.voxel_x param.voxel_z/param.interp];
        s = size(XguessSAVE1);
        centroid_orig = s/2.*[param.voxel_y param.voxel_x param.voxel_z];
        centroid_error = centroid_interp - centroid_orig;
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1);
        msg = sprintf('      original centroid (um) = [ %4.1f %4.1f %4.1f ]\n',centroid_orig(1),centroid_orig(2),centroid_orig(3));
        msg = sprintf('%sinterpolated centroid (um) = [ %4.1f %4.1f %4.1f ]\n',msg,centroid_interp(1),centroid_interp(2),centroid_interp(3));
        msg = sprintf('%s  centroid difference (um) = [ %4.1f %4.1f %4.1f ]\n',msg,centroid_error(1),centroid_error(2),centroid_error(3));
        text(0.2,0.8,msg,'FontSize',12,'Interpreter','none');
        str=sprintf('%s%s_intensity_projection_dim3_%s.png',param.savePath,param.timestamp,mystr);
        save_plot(f, str);
    end
else
    fprintf('File not found:\n%s\n\n',f);
    keyboard;
end 
end
