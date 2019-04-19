function GRID = updateGRID_flake(GRID)


%update ice cover and water grid
GRID.lake.ice.cT_domain = (GRID.general.K_grid(2:end) - GRID.general.K_grid(min([GRID.lake.water.cT_domain_ub GRID.lake.ice.cT_domain_ub GRID.lake.ice.cT_domain_ub])))<=GRID.lake.ice.z_ice & ...
    (GRID.general.K_grid(2:end) - GRID.general.K_grid(min([GRID.lake.water.cT_domain_ub GRID.lake.ice.cT_domain_ub])))> 0;

GRID.lake.ice.K_domain = GRID.lake.ice.cT_domain;
[GRID.lake.ice.cT_domain_lb GRID.lake.ice.cT_domain_ub] = LayerIndex(GRID.lake.ice.cT_domain);
[GRID.lake.ice.K_domain_lb GRID.lake.ice.K_domain_ub]   = LayerIndex(GRID.lake.ice.K_domain);

%update water body grid
GRID.lake.water.cT_domain = (~GRID.air.cT_domain &  ~GRID.snow.cT_domain & ~GRID.lake.ice.cT_domain & ~GRID.soil.cT_domain);
GRID.lake.water.K_domain  = GRID.lake.water.cT_domain;
[GRID.lake.water.cT_domain_lb GRID.lake.water.cT_domain_ub] = LayerIndex(GRID.lake.water.cT_domain);
[GRID.lake.water.K_domain_lb GRID.lake.water.K_domain_ub]   = LayerIndex(GRID.lake.water.K_domain);
