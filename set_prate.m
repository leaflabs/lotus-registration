function param = set_prate (MI, param)
% aim for 0.95 aceptance at T0
%param.T0 = -1/log(param.init_p);
%param.TC0 = round(log10(param.final_p) / log10(param.Trate))
%param.T0 = -1/log(0.95)/MI;
param.prate = 10^( ( log10(param.final_p)-log10(param.init_p) ) / param.TC0);
end