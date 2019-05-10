function [ inf_altitude, inf_cT_index, pfTable_altitude] = getInfiltrationAltitude( PARA, GRID, T)
% Function that gives the cT index of the bottom cell of the bucket.
inf_altitude = NaN;
inf_cT_index = NaN;
pfTable_altitude = NaN;

if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0; % Condition to work on infiltration
    T=T(GRID.soil.cT_domain);
    i=1;
    [~,i_max]=min(abs((PARA.location.altitude - GRID.general.K_grid(GRID.soil.K_domain))-PARA.soil.infiltration_limit_altitude));

    while  T(i)>0 && i<=i_max
        i=i+1;
    end
    inf_cT_index=i-1;
    inf_altitude = PARA.location.initial_altitude - PARA.location.shrinkage - GRID.general.K_grid(GRID.soil.cT_domain_ub+inf_cT_index);
    if i<i_max;
       pfTable_altitude=inf_altitude;
    end
end


end
