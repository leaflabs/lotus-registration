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
