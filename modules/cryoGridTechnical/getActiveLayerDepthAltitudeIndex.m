function [ ald_cT_index ] = getActiveLayerDepthAltitudeIndex( GRID, T)
% Function that gives the cT index of the bottom cell of the bucket.
    ald_cT_index = NaN;
    
    if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0; % Condition to work on infiltration
        T=T(GRID.soil.cT_domain);
        i=1;
        i_max=200;
        while  T(i)>0 && i<=i_max
            i=i+1;
        end
        ald_cT_index=GRID.soil.cT_domain_ub+i-2;
    end
end
