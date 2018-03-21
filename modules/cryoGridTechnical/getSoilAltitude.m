function soil_altitude = getSoilAltitude(PARA, GRID)

if ~isempty( GRID.lake.cT_domain_ub )
    soil_altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.lake.cT_domain_lb+1);
else
    soil_altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.soil.cT_domain_ub);
end