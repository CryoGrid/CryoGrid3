function [ ald_altitude ] = getActiveLayerDepthAltitude( PARA, GRID, T)
% Function that gives the cT index of the bottom cell of the bucket.
    ald_altitude = NaN;
    
    if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0; % Condition to work on infiltration
        T=T(GRID.soil.cT_domain);
        i=1;
        i_max=200;
        while  T(i)>0 && i<=i_max
            i=i+1;
        end
        ald_cT_index=i-1;
        ald_altitude = PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub+ald_cT_index);
    end
    

end
