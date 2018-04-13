function soil_altitude = getSoilAltitude(PARA, GRID)

%tsvd GRID.lake.cT_domain replaced by GRID.lake.water.cT_domain
if ~isempty( GRID.lake.water.cT_domain_ub )   %zzz rm what is soil_altitude needed for...? check....
    soil_altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.lake.water.cT_domain_lb+1);
else
    soil_altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.soil.cT_domain_ub);
end