    % at this point the thermal and hydrological state of the soil and snow is calculated
    % energy content soil domain in [J/m^2]
    % distinguished by sensible and latent part
    E_soil_sens_old = E_soil_sens;
    E_soil_lat_old = E_soil_lat;
    E_soil_old = E_soil;
    E_soil_sens = nansum(  ( PARA.constants.c_w .* lwc + PARA.constants.c_i .* (wc-lwc) + PARA.constants.c_m .* GRID.soil.cT_mineral + PARA.constants.c_o .* GRID.soil.cT_organic ) .* T(GRID.soil.cT_domain) ...
                       .* GRID.general.K_delta(GRID.soil.cT_domain)  );
    E_soil_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc .* GRID.general.K_delta(GRID.soil.cT_domain)  );
    E_soil = E_soil_sens + E_soil_lat;
    TEMPORARY.dE_soil_sens = TEMPORARY.dE_soil_sens + E_soil_sens - E_soil_sens_old;
    TEMPORARY.dE_soil_lat = TEMPORARY.dE_soil_lat + E_soil_lat - E_soil_lat_old;
    TEMPORARY.dE_soil = TEMPORARY.dE_soil + E_soil - E_soil_old;
    % energy content snow domain in [J/m^2]
    E_snow_sens_old = E_snow_sens;
    E_snow_lat_old = E_snow_lat;
    E_snow_old = E_snow;
    E_snow_sens = nansum( c_cTgrid(GRID.snow.cT_domain) .* T(GRID.snow.cT_domain) .* GRID.general.K_delta(GRID.snow.cT_domain) );
    E_snow_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.snow.cT_domain).* GRID.general.K_delta(GRID.snow.cT_domain) );
    E_snow = E_snow_sens + E_snow_lat;
    TEMPORARY.dE_snow_sens = TEMPORARY.dE_snow_sens + E_snow_sens - E_snow_sens_old;
    TEMPORARY.dE_snow_lat = TEMPORARY.dE_snow_lat + E_snow_lat - E_snow_lat_old;
    TEMPORARY.dE_snow = TEMPORARY.dE_snow + E_snow - E_snow_old;    
    % water content soil domain in [m]
    W_soil_old = W_soil;
    W_soil = nansum( wc .* GRID.general.K_delta(GRID.soil.cT_domain) );
    HYDRO.dW_soil = HYDRO.dW_soil + (W_soil - W_soil_old)*1000; % in [mm]    
    % water content snow domain in [m]
    W_snow_old = W_snow;
    W_snow = nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;
    HYDRO.dW_snow = HYDRO.dW_snow + (W_snow - W_snow_old)*1000; % in [mm]
  