function [  sediment_change_tot,...
            sediment_change_diff,...
            sediment_change_adv,...
            sediment_change_o,...
            sediment_change_m,...
            GRID                      ] = calculateLateralSedimentFluxes( PARA, GRID, PACKAGE_sedimentExchange_j, T, j )

sediment_change_tot = 0;
sediment_change_diff =0;
sediment_change_adv = 0;
sediment_change_o = 0;
sediment_change_m = 0;


if PARA.ensemble.distanceBetweenPoints(labindex,j)>0
    
    % conditions for erosion met in both realizations
    erosion_condition_index = isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0;
    erosion_condition_j = PACKAGE_sedimentExchange_j.erosion_condition;
    
    if erosion_condition_index && erosion_condition_j
        
        
        soil_surface_index = PARA.ensemble.soil_altitude(labindex);
        soil_surface_j = PARA.ensemble.soil_altitude(j);
        
        surface_index = PARA.ensemble.altitude(labindex);
        surface_j = PARA.ensemble.altitude(j);
        
        cT_mineral_j = PACKAGE_sedimentExchange_j.cT_mineral ;
        cT_organic_j = PACKAGE_sedimentExchange_j.cT_organic ;
        K_delta_j = PACKAGE_sedimentExchange_j.K_delta;
        area_j = PARA.ensemble.area(j);
        
        K_delta_index = GRID.general.K_delta(GRID.soil.cT_domain);
        area_index = PARA.ensemble.area(labindex);
        
        % set hillslope diffusivities and critical angle
        K_land = PARA.soil.hillslope_diffusivity_land;
        K_water = PARA.soil.hillslope_diffusivity_water;
        alpha_crit = PARA.soil.critical_hillslope_angle; % to be chosen depending on whether surface frozen or unfrozen, for now 45°
        
        % weighting of diffusive and advective transport
        w_diff = PARA.soil.weight_diffusion;
        w_adv = PARA.soil.weight_advection;
        
        % determine Keff for mixed interface
        phi_tot = abs( soil_surface_index - soil_surface_j );
        phi_land = max( soil_surface_index, soil_surface_j ) - min( surface_index, surface_j );
        phi_land = max( phi_land, 0 );
        phi_water = min( surface_index, surface_j ) - min( soil_surface_index, soil_surface_j);
        phi_water = max ( min( phi_water, phi_tot), 0 );
        assert( abs( phi_land + phi_water - phi_tot)<1e-9 , 'CryoGridLateralErosion - water/air interfaces do not match total interface' ) ;
        K_eff = phi_tot ./ ( phi_land./K_land + phi_water./K_water );       % based on assuming "series junction" of erosion domains --> reciprocal addition of "conductivities"
        
        D = PARA.ensemble.distanceBetweenPoints(labindex,j);
        L = PARA.ensemble.thermal_contact_length(labindex,j);
        
        % caclulate sediment fluxes due to diffusion
        q_sed_diff = w_diff .* K_eff .* (soil_surface_j - soil_surface_index) ./ D ;  % sediment flux in [m^3/ (m sec)]
        
        % calculate sediment fluxes due to advection
        grad = (soil_surface_j - soil_surface_index) ./ D;
        
        grad = sign(grad) .* min( abs(grad), 0.95.*tan(alpha_crit) ); % this is just a workaround to limit fluxes when they approach the critical angle
        
        
        assert( atan(grad)^2<alpha_crit^2, 'CryoGridLateralErosion - gradient exceeds critical slope angle' );
        q_sed_adv = w_adv .* K_eff .* grad .* atan(grad)^2 ./ ( alpha_crit^2 - atan(grad)^2 );
        
        % calculate total sediment volume transported within timestep
        V_sed_diff = abs( q_sed_diff .* L .* PARA.technical.syncTimeStep .* 24 .* 3600 );  % sediment volume in [m^3] (always positive)
        V_sed_adv = abs( q_sed_adv .* L .* PARA.technical.syncTimeStep .* 24 .* 3600 );  % sediment volume in [m^3] (always positive)
        
        % combine diffusive and advective fluxes
        q_sed = q_sed_diff + q_sed_adv;
        V_sed = V_sed_diff + V_sed_adv;
        
        % here, the total sediment flux should be limited somehow
        % ...
        % -------------
        
        if q_sed>0     % gaining sediment
            
            % calculate composition of organic and mineral
            V_sed_o = 0;
            V_sed_m = 0;
            V_sed_remaining = V_sed;
            k=1;
            while V_sed_remaining > 0
                if cT_organic_j(k)+cT_mineral_j(k)>0
                    V_temp = min( V_sed_remaining, K_delta_j(k) .* (cT_organic_j(k)+cT_mineral_j(k)) .* area_j );
                    V_sed_o = V_sed_o + cT_organic_j(k) ./ (cT_organic_j(k)+cT_mineral_j(k)) .* V_temp;
                    V_sed_m = V_sed_m + cT_mineral_j(k) ./ (cT_organic_j(k)+cT_mineral_j(k)) .* V_temp;
                    V_sed_remaining = V_sed_remaining - V_temp;
                end
                k=k+1;
            end
            
            sediment_change_tot = V_sed ./ area_index;
            sediment_change_diff = w_diff .* V_sed_diff ./ area_index;
            sediment_change_adv = w_adv .* V_sed_adv ./ area_index;
            sediment_change_o = V_sed_o ./ area_index;
            sediment_change_m = V_sed_m ./ area_index;
            
            assert( abs( V_sed_o + V_sed_m - V_sed )< 1e-12, 'lateral erosion - organic and mineral do not match total sed change' );
            assert( abs( V_sed_diff + V_sed_adv - V_sed )< 1e-12, 'lateral erosion - diffusive and advective do not match total sed change' );
            
        elseif q_sed<0     % losing sediment
            
            
            % calculate composition of organic and mineral
            V_sed_o = 0;
            V_sed_m = 0;
            V_sed_remaining = V_sed;
            k=1;
            while V_sed_remaining > 0
                if GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)>0
                    V_temp = min( V_sed_remaining, K_delta_index(k) .* (GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)) .* area_index );
                    V_sed_o = V_sed_o + GRID.soil.cT_organic(k) ./ (GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)) .* V_temp;
                    V_sed_m = V_sed_m + GRID.soil.cT_mineral(k) ./ (GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)) .* V_temp;
                    V_sed_remaining = V_sed_remaining - V_temp;
                end
                k=k+1;
            end
            
            
            sediment_change_tot = -V_sed ./ area_index;
            sediment_change_diff = -w_diff .* V_sed_diff ./ area_index;
            sediment_change_adv = -w_adv .* V_sed_adv ./ area_index;
            sediment_change_o = -V_sed_o ./ area_index;
            sediment_change_m = -V_sed_m ./ area_index;
            
            assert( abs( V_sed_o + V_sed_m - V_sed )< 1e-12, 'lateral erosion - organic and mineral do not match total sed change' );
            assert( abs( V_sed_diff + V_sed_adv - V_sed )< 1e-12, 'lateral erosion - diffusive and advective do not match total sed change' );
            
        else % same altitude
            ...
        end
    else % no erosion condition
        ...
    end
else % not connected
    ...
end