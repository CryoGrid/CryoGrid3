function GRID=makeGrids(PARA)

    GRID.snow.snowCellSize=PARA.technical.SWEperCell/(PARA.snow.rho_snow/PARA.constants.rho_w);
    GRID.snow.snowGrid=[-1.*(PARA.technical.maxSWE./(PARA.technical.SWEperCell)+2).*GRID.snow.snowCellSize:GRID.snow.snowCellSize:-GRID.snow.snowCellSize]';
    GRID.soil.soilGrid=PARA.technical.subsurfaceGrid;

    K_grid =[GRID.snow.snowGrid;  GRID.soil.soilGrid]; %grid on which the conductivty information lives (edges of grid cells)
    cT_grid=(K_grid(1:end-1)+K_grid(2:end))/2; %grid on which heat capacity and temperature information lives (midpoints of grid cells)
    cT_delta=(-cT_grid(1:end-1,1)+cT_grid(2:end,1));
    K_delta=(-K_grid(1:end-1,1)+K_grid(2:end,1));

    GRID.general.cT_grid=cT_grid;
    GRID.general.K_grid=K_grid;
    GRID.general.cT_delta=cT_delta;
    GRID.general.K_delta=K_delta;

    %set air grid
    GRID.air.cT_domain= false(size(GRID.general.cT_grid));
    GRID.air.cT_domain(1:length(GRID.snow.snowGrid))=1;
    GRID.air.K_domain= false(size(GRID.general.K_grid));
    GRID.air.K_domain(1:length(GRID.snow.snowGrid))=1;
    [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
    [GRID.air.K_domain_lb, GRID.air.K_domain_ub]   = LayerIndex(GRID.air.K_domain);

    %set soil grid
    GRID.soil.cT_domain= false(size(GRID.general.cT_grid));
    GRID.soil.cT_domain(GRID.general.cT_grid>0)=1;
    [GRID.soil.cT_domain_lb, GRID.soil.cT_domain_ub] = LayerIndex(GRID.soil.cT_domain);
    GRID.soil.K_domain= false(size(GRID.general.K_grid));
    GRID.soil.K_domain(GRID.soil.cT_domain_ub:end)=1;
    [GRID.soil.K_domain_lb, GRID.soil.K_domain_ub]   = LayerIndex(GRID.soil.K_domain);

    %set snow grid
    GRID.snow.cT_domain=false(size(GRID.general.cT_grid));
    GRID.snow.K_domain=false(size(GRID.general.K_grid));
    [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
    [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub]   = LayerIndex(GRID.snow.K_domain);

end