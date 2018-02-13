function [ ald_altitude ] = getActiveLayerDepthAltitude( PARA, GRID, T,index)
% Function that gives the cT index of the bottom cell of the bucket.
ald_altitude = NaN;




if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0; % Condition to work on infiltration
    T=T(GRID.soil.cT_domain);
    i=1;
    
    if numlabs < 2;
        [~,i_max]=min(abs((PARA.location.altitude(index)-GRID.soil.soilGrid)-PARA.soil.alt_infiltration_limit));
    else
        [~,i_max]=min(abs((PARA.ensemble.altitude(index)-GRID.soil.soilGrid)-PARA.ensemble.alt_infiltration_limit));
    end
    
    % i_max=length(T)-1; %200;
    while  T(i)>0 && i<=i_max
        i=i+1;
    end
    ald_cT_index=i-1;
    ald_altitude = PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub+ald_cT_index);
end


end
