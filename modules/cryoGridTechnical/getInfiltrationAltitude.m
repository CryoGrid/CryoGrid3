function [ inf_altitude, inf_cT_index] = getInfiltrationAltitude( PARA, GRID, T)
% Function that gives the cT index of the bottom cell of the bucket.
inf_altitude = NaN;
inf_cT_index = NaN;

if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0; % Condition to work on infiltration
    T=T(GRID.soil.cT_domain);
    i=1;
    [~,i_max]=min(abs((PARA.location.altitude - GRID.soil.soilGrid)-PARA.soil.infiltration_limit_altitude));

<<<<<<< HEAD
    
    % i_max=length(T)-1; %200;
=======
>>>>>>> origin/xice_mpi_polygon_TC
    while  T(i)>0 && i<=i_max
        i=i+1;
    end
    inf_cT_index=i-1;
    inf_altitude = PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub+inf_cT_index);
end


end
