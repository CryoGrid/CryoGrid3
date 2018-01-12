function BALANCE = updateBALANCE(T, wc, c_cTgrid, lwc_cTgrid, BALANCE, GRID, PARA)

    % at this point the thermal and hydrological state of the soil and snow is calculated
    % energy content soil domain in [J/m^2]
    % distinguished by sensible and latent part
    E_soil_sens_old = BALANCE.energy.E_soil_sens;
    E_soil_lat_old = BALANCE.energy.E_soil_lat;
    E_soil_old = BALANCE.energy.E_soil;
    BALANCE.energy.E_soil_sens = nansum(  ( PARA.constants.c_w .* lwc_cTgrid(GRID.soil.cT_domain) + PARA.constants.c_i .* (wc-lwc_cTgrid(GRID.soil.cT_domain)) + PARA.constants.c_m .* GRID.soil.cT_mineral + PARA.constants.c_o .* GRID.soil.cT_organic ) .* T(GRID.soil.cT_domain) ...
                       .* GRID.general.K_delta(GRID.soil.cT_domain)  );
    BALANCE.energy.E_soil_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.soil.cT_domain) .* GRID.general.K_delta(GRID.soil.cT_domain)  );
    BALANCE.energy.E_soil = BALANCE.energy.E_soil_sens + BALANCE.energy.E_soil_lat;
    BALANCE.energy.dE_soil_sens = BALANCE.energy.dE_soil_sens + BALANCE.energy.E_soil_sens - E_soil_sens_old;
    BALANCE.energy.dE_soil_lat = BALANCE.energy.dE_soil_lat + BALANCE.energy.E_soil_lat - E_soil_lat_old;
    BALANCE.energy.dE_soil = BALANCE.energy.dE_soil + BALANCE.energy.E_soil - E_soil_old;
    % energy content snow domain in [J/m^2]
    E_snow_sens_old = BALANCE.energy.E_snow_sens;
    E_snow_lat_old = BALANCE.energy.E_snow_lat;
    E_snow_old = BALANCE.energy.E_snow;
    BALANCE.energy.E_snow_sens = nansum( c_cTgrid(GRID.snow.cT_domain) .* T(GRID.snow.cT_domain) .* GRID.general.K_delta(GRID.snow.cT_domain) );
    BALANCE.energy.E_snow_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.snow.cT_domain).* GRID.general.K_delta(GRID.snow.cT_domain) );
    BALANCE.energy.E_snow = BALANCE.energy.E_snow_sens + BALANCE.energy.E_snow_lat;
    BALANCE.energy.dE_snow_sens = BALANCE.energy.dE_snow_sens + BALANCE.energy.E_snow_sens - E_snow_sens_old;
    BALANCE.energy.dE_snow_lat = BALANCE.energy.dE_snow_lat + BALANCE.energy.E_snow_lat - E_snow_lat_old;
    BALANCE.energy.dE_snow = BALANCE.energy.dE_snow + BALANCE.energy.E_snow - E_snow_old;    
    
    % water content soil domain in [m]
    W_soil_old = BALANCE.water.W_soil;
    BALANCE.water.W_soil = nansum( wc .* GRID.general.K_delta(GRID.soil.cT_domain) );
    BALANCE.water.dW_soil = BALANCE.water.dW_soil + (BALANCE.water.W_soil - W_soil_old)*1000; % in [mm]    
    % water content snow domain in [m]
    W_snow_old = BALANCE.water.W_snow;
    BALANCE.water.W_snow = nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;
    BALANCE.water.dW_snow = BALANCE.water.dW_snow + (BALANCE.water.W_snow - W_snow_old)*1000; % in [mm]
    
end