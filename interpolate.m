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
    debug_count = 0;
    for i=1:N
        if i < boundary
            if i>1
                out(:,:,i-1) = LFM(:,:,1);
            end
        elseif i > (N-(boundary-1))
            out(:,:,i-1) = LFM(:,:,end);
        else
            v = a*A+b*B;
            fprintf('count = %d, a = %1.3f, b = %1.3f, a+b = %1.3f, j = %d, j+1 = %d, A = %d, B = %d, v = %d\n',debug_count,a,b,a+b,j,j+1,uint32(sum(sum(A))),uint32(sum(sum(B))),uint32(sum(sum(v))));
            if debug_count==306
                keyboard
            end
            debug_count = debug_count + 1;
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
