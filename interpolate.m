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