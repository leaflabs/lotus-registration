function bead_center = axialPSF (config, handle, axialVals, crop_index, dim, b, z, myflag)
myStr = {};
threshold = config.threshold;
N = length(axialVals)-1;
axialPosInt0 = 0;
axialPosInt1 = 0;
riseFlag = false;
fallFlag = false;
zvals = [z-5:z+5];
% JPK
offset = crop_index(3) * config.zspacing;
%
% for each axialVals = number of columns in cropped image
for j=1:length(axialVals)-1
    axialVals(j);
    if j<config.axialBuf | N-j+1<config.axialBuf
        %j
        continue;
    end
    if axialVals(j) <= threshold && axialVals(j+1) > threshold
        if myflag
            keyboard
        end
        if ~riseFlag
            riseFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max axial spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            axialPosInt0 = -1;
            %             % calculate bead center
            %             bead_center = [-1 -1];
            %             thisStr = ['Setting bead axial center and full-width half-max axial spread to %f and %f.', ...
            %                 bead_center(1),bead_center(2)];
            %             disp(thisStr);
            %             return;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        %
        % crop_ind = [minrow maxrow mincol maxcol]
        axialPos0 = offset + j*config.zspacing;
        %axialPos0 = allLabels(zvals(j));
        axialPos1 = axialPos0 + config.zspacing;
        scale = ( threshold - axialVals(j) ) / ( axialVals(j+1) - axialVals(j) );
        axialPosInt0 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        %axialPosInt0 = (j+scale)*pixel;
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(axialPosInt0))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    elseif axialVals(j) >= threshold && axialVals(j+1) < threshold
        if myflag
            keyboard
        end
        if ~fallFlag
            fallFlag = true;
        else
            thisStr = ['Warning! Too many falling threshold crossings! Full-width half-max axial spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            axialPosInt0 = -1;
            %             % calculate bead center
            %             bead_center = [-1 -1];
            %             thisStr = ['Setting bead axial center and full-width half-max axial spread to %f and %f.', ...
            %                 bead_center(1),bead_center(2)];
            %             disp(thisStr);
            %             return;
        end
        if axialPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        axialPos0 = offset + j*config.zspacing;
        %axialPos0 = allLabels(zvals(j));
        axialPos1 = axialPos0 + config.zspacing;
        %axialPos0 =  allLabels(zvals(j));
        %axialPos1 = allLabels(zvals(j+1));
        scale = ( axialVals(j) - threshold) / ( axialVals(j) - axialVals(j+1) );
        axialPosInt1 = axialPos0 + (axialPos1-axialPos0)*scale; %Equivalently...
        %axialPosInt1 = (j+scale)*pixel;
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(axialPosInt1))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    end
end
if myflag
    keyboard
end
if axialPosInt0 < 0
    thisStr = ['Warning! Full-width half-max axial spread not calculated.'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    myStr = myStr(3:end);
    % calculate bead center
    bead_center = [-1 -1];
    %thisStr = ['Setting bead axial center and full-width half-max axial spread to %f and %f.', ...
    %    bead_center(1),bead_center(2)];
    %disp(thisStr);
else
    if myflag
        keyboard
    end
    fwhmA = axialPosInt1-axialPosInt0;
    thisStr = ['Full-width half-max axial spread = ' num2str(round(fwhmA)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = [0.5*(axialPosInt1+axialPosInt0) fwhmA];
end
figure(handle);
x = xlim;
y = ylim;
hold on;
plot(xlim,(y(1)+config.latBuf-1)*[1 1],'r-');
plot(xlim,(y(2)-config.latBuf+1)*[1 1],'r-');
hold off;
subplot(2,2,3);
text(5.0,0.5,myStr,'FontSize',6);
%% duplicate
thisdim = get(gca,'Position');
set(gca, 'Position', [thisdim(1) thisdim(2) dim(3) thisdim(4)]);
set(gca,'XTickLabel',{});
xlim([1 b(2)]);
%
%
x = xlim;
y = ylim;
hold on;
plot((x(1)+config.axialBuf-1)*[1 1],ylim,'r-');
plot((x(2)-config.axialBuf+1)*[1 1],ylim,'r-');
hold off;
end