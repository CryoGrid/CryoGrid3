function surface_altitude = getSurfaceAltitude(PARA, GRID)

surface_altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.air.cT_domain_lb+1);