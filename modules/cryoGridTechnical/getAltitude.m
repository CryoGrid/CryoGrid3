function altitude = getAltitude( PARA, GRID )
    
    altitude = PARA.location.initial_altitude - PARA.location.shrinkage - GRID.general.K_grid(GRID.soil.cT_domain_ub);

end