function water_fluxes = calculateLateralSurfaceSubsurfaceWaterDarcyFluxes(T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes)
% Function that calculate the water fluxes between the 2 workers.
% PARA.ensemble should contain the elevations, the water table values, the
% bucketbottom the weight and conductivities of each worcker.
% GRID will contain the soil water content.
% T is the temperature vector of the current worker.
% packageWorkerj has to be a bundle as defined considering necessary inputs
index = labindex;

% Check whether the considered workers (index, j) are hydrologically connceted
if PARA.ensemble.adjacency_water(index,j)>0
    % conditions for water infiltration / exchanges
    infiltration_condition_index = isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0;
    infiltration_condition_j = PACKAGE_waterExchange_j.infiltration_condition;
    
    if infiltration_condition_index && infiltration_condition_j
        
        wt_index = PARA.ensemble.water_table_altitude(index);
        wt_j     = PACKAGE_waterExchange_j.water_table_altitude;
        inf_index = PARA.ensemble.infiltration_altitude(index);
        inf_j     = PACKAGE_waterExchange_j.infiltration_altitude;
        soil_index = PARA.ensemble.soil_altitude(index);
        soil_j = PACKAGE_waterExchange_j.soil_altitude;

        % Decipher between cases
        [waterpotWindex, hasWater_index] = nanmax([wt_index, inf_index ] );
        [waterpotWj, hasWater_j]         = nanmax([wt_j,     inf_j     ] );
        
        waterpotW_max = max( [ waterpotWindex, waterpotWj ] );
        inf_max = max( [ inf_index, inf_j ] );
        soil_max = max( [ soil_index, soil_j ] );        
        
        % determine surface and subsurface flux domains
        H_surf = max( [ 0, waterpotW_max - max( [ soil_max, inf_max ] ) ] );
        H_subs = max( [ 0, min( [ waterpotW_max - inf_max, soil_max - inf_max ] ) ] );
        
        % check consistency
        assert( abs( H_surf+H_subs - max( [ 0, waterpotW_max - inf_max ] ) ) < 1e-9 , 'contact heights do not match' );
       
        if (waterpotWj > waterpotWindex && hasWater_j==1)% Current worker is gaining water
            % gradient of water table which is used for both surface and subsurface fluxes
            gradW = (waterpotWj - waterpotWindex) ./ PARA.ensemble.hydraulicDistance(index,j) ;
       
            % calculate fluxes through surface/subsurface domain cross sections in [m^3/s]
            FluxSurface = PARA.ensemble.hydraulic_conductivity_surf(index,j) .* gradW .* H_surf .* PARA.ensemble.hydraulic_contact_length(index,j);
                
            % calcuate subsurface fluxes according to Darcy law
            FluxSubsurface = PARA.ensemble.hydraulic_conductivity_subs(index,j) .* gradW .* H_subs .* PARA.ensemble.hydraulic_contact_length(index,j);

            % calcualte height changes in entire sync time interval
            water_fluxes(index,j) = (FluxSurface+FluxSubsurface) * PARA.technical.syncTimeStep * 3600 * 24 / PARA.ensemble.area(index);        
            %water_fluxes(j,index) = -1 * (DarcyFluxSurface+DarcyFluxSubsurface) * PARA.technical.syncTimeStep * 3600 * 24 / PARA.ensemble.area(j);
            
        elseif (waterpotWindex > waterpotWj && hasWater_index==1) % Current worker is loosing water
            % gradient of water table which is used for both surface and subsurface fluxes
            gradW = (waterpotWindex - waterpotWj) ./ PARA.ensemble.hydraulicDistance(index,j) ;
            
            % calculate fluxes through surface/subsurface domain cross sections in [m^3/s]
            FluxSurface = PARA.ensemble.hydraulic_conductivity_surf(index,j) .* gradW .* H_surf .* PARA.ensemble.hydraulic_contact_length(index,j);
            
            % calcuate subsurface fluxes according to Darcy law
            FluxSubsurface = PARA.ensemble.hydraulic_conductivity_subs(index,j) .* gradW .* H_subs .* PARA.ensemble.hydraulic_contact_length(index,j);

            water_fluxes(index,j) = -1 * (FluxSurface+FluxSubsurface) * PARA.technical.syncTimeStep * 3600 * 24 / PARA.ensemble.area(index);        
            %water_fluxes(j,index) = (DarcyFluxSurface+DarcyFluxSubsurface) * PARA.technical.syncTimeStep * 3600 * 24 / PARA.ensemble.area(j);       
                       
        else % same water table, no lateralFlux
            water_fluxes(j,index)=0;
            water_fluxes(index,j)=0;
        end
        
    else % No possible movement of water, because of snow or suficial freezing
        water_fluxes(j,index)=0;
        water_fluxes(index,j)=0;
    end
    
else % No connection between workers
    water_fluxes(j,index) = NaN;
    water_fluxes(index,j) = NaN;
end

end
