function [T, snow_i, snow_w]=refreeze(T, snow_i, snow_w, snow_a, c_temp, PARA)

    L=PARA.constants.L_sl.*PARA.constants.rho_w;%3.34e8;

    delta_SWE = min([snow_w,  -1.*double(T<=0).*T.*c_temp.*(snow_i+snow_w+snow_a)./L]'); 
    delta_SWE=delta_SWE';

    snow_i=snow_i + delta_SWE;
    snow_w=snow_w - delta_SWE;

    T=T + delta_SWE.*L ./ ((snow_i+snow_w+snow_a).*c_temp);

    T(isnan(T(:,1)),1)=0;


end