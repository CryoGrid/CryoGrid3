function waterTable = getWaterTableAltitude(T, wc, GRID, PARA)

    T=T(GRID.soil.cT_domain);
    K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m

    if ~isempty(GRID.snow.cT_domain_ub) || T(1)<=0
        waterTable=NaN;
    else
        waterTable = PARA.location.altitude;
        water = wc;
        porosity=(1. - GRID.soil.cT_mineral - GRID.soil.cT_organic);
        i=1;
        while water(i)<porosity(i) && T(i)>0
            waterTable = waterTable - K_delta(i);
            i=i+1;
        end
        
        if T(i)<=0;
            waterTable=NaN;
            disp('Dry to permafrost - waterTable=NaN \n')% Loop was stopped by the temperature condition so no water table
        end
    end

end
