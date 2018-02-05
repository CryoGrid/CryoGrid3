function surface_altitude = getSurfaceAltitude(PARA, GRID)

surface_altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.air.cT_domain_lb+1);

%if ~isempty(GRID.snow.cT_domain_ub)
%    surface_altitude = PARA.location.altitude + sum( GRID.general.K_delta(GRID.snow.cT_domain) );
%else
%    surface_altitude = PARA.location.altitude;
%end
    
