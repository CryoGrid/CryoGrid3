% energy content soil domain in [J/m^2]
% distinguished by sensible and latent part
E_soil_sens = nansum(  ( PARA.constants.c_w .* lwc + PARA.constants.c_i .* (wc-lwc) + PARA.constants.c_m .* GRID.soil.cT_mineral + PARA.constants.c_o .* GRID.soil.cT_organic ) .* T(GRID.soil.cT_domain) ...
                   .* GRID.general.K_delta(GRID.soil.cT_domain)  );
E_soil_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc .* GRID.general.K_delta(GRID.soil.cT_domain)  );
E_soil = E_soil_sens + E_soil_lat;
% energy content snow domain in [J/m^2]
E_snow_sens = nansum( c_cTgrid(GRID.snow.cT_domain) .* T(GRID.snow.cT_domain) .* GRID.general.K_delta(GRID.snow.cT_domain) );
E_snow_lat = nansum( PARA.constants.rho_w .* PARA.constants.L_sl .* lwc_cTgrid(GRID.snow.cT_domain).* GRID.general.K_delta(GRID.snow.cT_domain) );
E_snow = E_snow_sens + E_snow_lat;
% water content soil domain in [m]
W_soil = nansum( wc .* GRID.general.K_delta(GRID.soil.cT_domain) );
% water content snow domain in [m]
W_snow = nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;