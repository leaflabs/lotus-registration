function bead_center = latPSF (config, handle, lateralVals)
myStr = {};
threshold = config.threshold;
N = length(lateralVals)-1;
latPosInt0 = 0;
latPosInt1 = 0;
riseFlag = false;
fallFlag = false;
for j=1:length(lateralVals)-1
    if j<config.latBuf | N-j+1<config.latBuf
        j
        continue;
    end
    if lateralVals(j) <= threshold && lateralVals(j+1) > threshold
        if ~riseFlag
            riseFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max lat spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            latPosInt0 = -1;
            % calculate bead center
            %             bead_center = [-1 -1];
            %             thisStr = ['Setting bead lat center and full-width half-max lat spread to %f and %f.', ...
            %                 bead_center(1),bead_center(2)];
            %             disp(thisStr);
            %             return;
        end
        if latPosInt0 < 0
            break;
        end
        % found a rising threshold crossing
        latPos0 =  j*config.pixel;
        latPos1 = (j+1)*config.pixel;
        scale = ( threshold - lateralVals(j) ) / ( lateralVals(j+1) - lateralVals(j) );
        %latPosInt = latPos0 + (latPos1-latPos0)*scale; Equivalently...
        latPosInt0 = (j+scale)*config.pixel;
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(latPosInt0))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    elseif lateralVals(j) >= threshold && lateralVals(j+1) < threshold
        if ~fallFlag
            fallFlag = true;
        else
            thisStr = ['Warning! Too many rising threshold crossings! Full-width half-max lat spread not calculated.'];
            disp(thisStr);
            myStr = [myStr;{thisStr}];
            latPosInt0 = -1;
            %             % calculate bead center
            %             bead_center = [-1 -1];
            %             thisStr = ['Setting bead lat center and full-width half-max lat spread to %f and %f.', ...
            %                 bead_center(1),bead_center(2)];
            %             disp(thisStr);
            %             return;
        end
        if latPosInt0 < 0
            break;
        end
        % found a falling threshold crossing
        latPos0 =  j*config.pixel;
        latPos1 = (j+1)*config.pixel;
        scale = ( lateralVals(j) - threshold) / ( lateralVals(j) - lateralVals(j+1) );
        %latPosInt = latPos0 + (latPos1-latPos0)*scale; Equivalently...
        latPosInt1 = (j+scale)*config.pixel;
        thisStr = ['Interpolated threshold (' sprintf('%g',threshold) ') = ' num2str(round(latPosInt1))];
        disp(thisStr);
        myStr = [myStr;{thisStr}];
    end
end
if latPosInt0 < 0
    thisStr = ['Warning! Full-width half-max lat spread not calculated.'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    myStr = myStr(3:end);
    % calculate bead center
    bead_center = [-1 -1];
    %     thisStr = ['Setting bead lat center and full-width half-max lat spread to %f and %f.', ...
    %         bead_center(1),bead_center(2)];
    %     disp(thisStr);
else
    fwhmL = latPosInt1-latPosInt0;
    thisStr = ['Full-width half-max lateral spread = ' num2str(round(fwhmL)) ' um'];
    disp(thisStr);
    myStr = [myStr;{thisStr}];
    % calculate bead center
    bead_center = [0.5*(latPosInt1+latPosInt0) fwhmL];
end

figure(handle);
subplot(2,2,2);
text(0.1,15.0,myStr,'FontSize',6);
%

end