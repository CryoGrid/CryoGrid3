function BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA)
    
    % ENERGY balance
    % energy content soil domain in [J/m^2] distinguished by sensible and latent part
    BALANCE.energy.E_soil_sens = nansum(  ( PARA.constants.c_w .* lwc_cTgrid(GRID.soil.cT_domain) + PARA.constants.c_i .* (wc-lwc_cTgrid(GRID.soil.cT_domain)) + PARA.constants.c_m .* GRID.soil.cT_mineral + PARA.constants.c_o .* GRID.soil.cT_organic ) .* T(GRID.soil.cT_domain) ...
                       .* GRID.general.K_delta(GRID.soil.cT_domain)  );
    BALANCE.energy.E_soil_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.soil.cT_domain) .* GRID.general.K_delta(GRID.soil.cT_domain)  );
    BALANCE.energy.E_soil = BALANCE.energy.E_soil_sens + BALANCE.energy.E_soil_lat;
    % energy content snow domain in [J/m^2] distinguished by sensible and latent part
    BALANCE.energy.E_snow_sens = nansum( c_cTgrid(GRID.snow.cT_domain) .* T(GRID.snow.cT_domain) .* GRID.general.K_delta(GRID.snow.cT_domain) );
    BALANCE.energy.E_snow_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.snow.cT_domain).* GRID.general.K_delta(GRID.snow.cT_domain) );
    BALANCE.energy.E_snow = BALANCE.energy.E_snow_sens + BALANCE.energy.E_snow_lat;
    % accumulated changes per output timestep
    BALANCE.energy.dE_soil_sens = 0;
    BALANCE.energy.dE_soil_lat = 0;
    BALANCE.energy.dE_soil = 0;  
    BALANCE.energy.dE_snow_sens = 0;
    BALANCE.energy.dE_snow_lat = 0;
    BALANCE.energy.dE_snow = 0;
    
    BALANCE.energy.Q_lateral = zeros( length(GRID.general.cT_grid) , 1 );
    
    
    % WATER balance
    % water content soil domain in [m]
    BALANCE.water.W_soil = nansum( wc .* GRID.general.K_delta(GRID.soil.cT_domain) );
    % water content snow domain in [m]
    BALANCE.water.W_snow = nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;
    % accumulated changes per output timestep
    % storage
    BALANCE.water.dW_soil = 0;
    BALANCE.water.dW_snow = 0;
    % precipitation
    BALANCE.water.dp_rain=0;
    BALANCE.water.dp_snow=0; % SWE
    % evapotranspiration and sublimation
    BALANCE.water.de=0;
    BALANCE.water.ds=0;
    % runoff
    BALANCE.water.dr_surface=0;
    BALANCE.water.dr_external=0;
    BALANCE.water.dr_snowmelt=0;
    BALANCE.water.dr_excessSnow=0;
    BALANCE.water.dr_lateralSnow=0;
    BALANCE.water.dr_rain=0;  % this is only rain on frozen ground
    BALANCE.water.dr_lateralWater=0;
    BALANCE.water.dr_DarcyReservoir=0; % When worker is connected to a Darcy_reservoir as a boundary condition
    BALANCE.water.dr_lateralExcess=0; % excess water when applying lateral fluxes
    % mismatch
    BALANCE.water.dm_lacking=0;
    
    if PARA.modules.lateral
        BALANCE.water.lateral_water_fluxes_vector=zeros(1,numlabs);
        BALANCE.water.lateral_water_fluxes_matrix=zeros(numlabs,numlabs);
        BALANCE.water.lateral_water_fluxes_boundary = 0;
        %BALANCE.water.dr_water_fluxes_out=zeros(numlabs,numlabs);
    end
    
    
    % MASS (SEDIMENT)
    BALANCE.mass.M_organic = nansum( GRID.general.K_delta(GRID.soil.cT_domain) .* GRID.soil.cT_organic ) + GRID.soil.residualOrganic;
    BALANCE.mass.M_mineral = nansum( GRID.general.K_delta(GRID.soil.cT_domain) .* GRID.soil.cT_mineral ) + GRID.soil.residualMineral;
    BALANCE.mass.dM_organic = 0;
    BALANCE.mass.dM_mineral = 0;
end