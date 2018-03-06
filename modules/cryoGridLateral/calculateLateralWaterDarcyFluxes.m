function water_fluxes = calculateLateralWaterDarcyFluxes(T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes)
% Function that calculate the water fluxes between the 2 workers.
% PARA.ensemble should contain the elevations, the water table values, the
% bucketbottom the weight and conductivities of each worcker.
% GRID will contain the soil water content.
% T is the temperature vector of the current worker.
% packageWorkerj has to be a bundle as defined considering necessary inputs

index = labindex;

% Check whether the considered workers (index, j) are hydrologically connceted

if PARA.ensemble.hydraulic_contact_length(index,j)>0
    % conditions for water infiltration / exchanges
    infiltration_condition_index = isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0;
    infiltration_condition_j = PACKAGE_waterExchange_j.infiltration_condition;
    
    if infiltration_condition_index && infiltration_condition_j
        
        wt_index = PARA.ensemble.water_table_altitude(index);
        wt_j     = PACKAGE_waterExchange_j.water_table_altitude;
        ald_index = PARA.ensemble.active_layer_depth_altitude(index);
        ald_j     = PACKAGE_waterExchange_j.active_layer_depth_altitude;
        

        % Decipher between cases
        [waterpotWindex, hasWater_index] = nanmax([wt_index, ald_index ] );
        [waterpotWj, hasWater_j]         = nanmax([wt_j,     ald_j     ] );
        
     
        if (waterpotWj > waterpotWindex && hasWater_j==1)% Current worker is gaining water

            % Calculate the maximum exchanged water volume
            DeltaH = waterpotWj - waterpotWindex;
            contact_height = waterpotWj - nanmax( [ waterpotWindex, ald_j ]);
            Distance=sqrt(DeltaH^2 + PARA.ensemble.hydraulicDistance(j,index)^2);
            section=contact_height .* PARA.ensemble.hydraulic_contact_length(j,index);
            DarcyFlux=PARA.ensemble.hydraulic_conductivity(j,index) * (DeltaH/Distance) * section; % in m3/sec
            % Attribute water height changes
            water_fluxes(j,index) = -1 * (DarcyFlux * PARA.technical.syncTimeStep * 3600 * 24 / PARA.ensemble.area(j)); % syncTimeStep in days, Darcy flux in m3/sec
            water_fluxes(index,j) = DarcyFlux * PARA.technical.syncTimeStep * 3600 * 24 / PARA.ensemble.area(index); % syncTimeStep in days, Darcy flux in m3/sec
           

        elseif (waterpotWindex > waterpotWj  && hasWater_index==1) % Current worker is loosing water

            % Calculate maximum of the exchange water volume
            DeltaH= waterpotWindex - waterpotWj;
            contact_height =waterpotWindex - nanmax( [ waterpotWj, ald_index ] );
            Distance=sqrt(DeltaH^2 + PARA.ensemble.hydraulicDistance(j,index)^2);
            section=contact_height .* PARA.ensemble.hydraulic_contact_length(j,index);
            DarcyFlux=PARA.ensemble.hydraulic_conductivity(j,index) * (DeltaH/Distance) * section; % in m3/sec
            % Attribute water height changes
            water_fluxes(j,index) = DarcyFlux * PARA.technical.syncTimeStep *24 *3600 / PARA.ensemble.area(j); % syncTimeStep in days, Darcy flux in m3/sec
            water_fluxes(index,j) = -1 * (DarcyFlux * PARA.technical.syncTimeStep *24 *3600 / PARA.ensemble.area(index)); % syncTimeStep in days, Darcy flux in m3/sec
            

        else % same water table, no lateralFlux
            water_fluxes(j,index)=0;
            water_fluxes(index,j)=0;
        end
        
    else % No possible movement of water, because of snow or suficial freezing
        water_fluxes(j,index)=0;
        water_fluxes(index,0)=0;
    end
    
else % No connection between workers
    water_fluxes(j,index) = NaN;
    water_fluxes(index,j) = NaN;
end

end
