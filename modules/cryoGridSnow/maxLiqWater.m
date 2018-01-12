function maxLiqWater=maxLiqWater(T, snow_i, snow_w, snow_a, poreSpace, c_temp, PARA)

    L=PARA.constants.L_sl.*PARA.constants.rho_w; %3.34e8;

    waterHoldingPot=0.05.* (poreSpace.*snow_i)./(1-poreSpace);   %in m; 5% of the pore space can be filled by water 

    maxLiqWater = (waterHoldingPot - snow_w - T.*c_temp.*(snow_i+snow_w+snow_a)./L);  %in m
    maxLiqWater = min([snow_a maxLiqWater]');
    maxLiqWater = maxLiqWater'; % negative if snow_w>waterHoldingPot

end