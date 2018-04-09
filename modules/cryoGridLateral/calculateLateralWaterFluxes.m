function water_flux_j = calculateLateralWaterFluxes(T, PACKAGE_waterExchange_j, GRID, PARA, j)
% Function that calculate the water fluxes between the 2 workers.
% PARA.ensemble should contain the elevations, the water table values, the
% bucketbottom the weight and conductivities of each worcker.
% GRID will contain the soil water content.
% T is the temperature vector of the current worker.
% packageWorkerj has to be a bundle as defined considering necessary inputs

index = labindex;

water_flux_j = NaN;


% Check whether the considered workers (index, j) are
% hydrologically connceted

if PARA.ensemble.hydraulic_contact_length(index,j)>0
    % conditions for water infiltration / exchanges
    infiltration_condition_index = isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0;
    infiltration_condition_j = PACKAGE_waterExchange_j.infiltration_condition;
    
    if infiltration_condition_index && infiltration_condition_j
        
        wt_index = PARA.ensemble.water_table_altitude(index);
        wt_j     = PACKAGE_waterExchange_j.water_table_altitude;
        inf_index = PARA.ensemble.infiltration_altitude(index);
        inf_j     = PACKAGE_waterExchange_j.infiltration_altitude;
        

        % Decipher between cases
        [waterpotWindex, hasWater_index] = nanmax([wt_index, inf_index ] );
        [waterpotWj, hasWater_j]         = nanmax([wt_j,     inf_j     ] );
        
     
        if (waterpotWj > waterpotWindex && hasWater_j==1)% Current worker is gaining water

            % Calculate the maximum exchanged water volume
            DeltaH = waterpotWj - waterpotWindex;
            contact_height = waterpotWj - nanmax( [ waterpotWindex, inf_j ]);
            Distance=sqrt(DeltaH^2 + PARA.ensemble.hydraulicDistance(j,index)^2);
            %section=(PARA.ensemble.waterTable(j)-soilGrid_workerj(PARA.ensemble.bottomBucketSoilcTIndex(j)+1)) * PARA.ensemble.contactLength(j,index);
            section=contact_height .* PARA.ensemble.hydraulic_contact_length(j,index);
            DarcyFlux=PARA.ensemble.hydraulic_conductivity(j,index) * (DeltaH/Distance) * section / PARA.ensemble.area(index);
            %MaxWaterHeight=DarcyFlux * PARA.technical.syncTimeStep ; % Jan : weight here corresponds to an actual area in m^2
            
%             % Calculate available water volume that can be lost
%             availableWater_li=Wc_workerj > PARA.soil.fieldCapacity; % Find where there is water above field capacity
%             cellresidual_water=0;
%             
%             if (altitude_j-soilGrid_workerj(PARA.ensemble.bottomBucketSoilcTIndex(j)+1))> waterpotWindex; % Bucket bottom of the worker loosing water is higher than water table of the worker gaining water
%                 
%                 availableWater_li(PARA.ensemble.bottomBucketSoilcTIndex(j)+1:end)=0; % Stop accounting for water after the bottom of the bucket
%                 
%             else % Bucket bottom of the worker loosing water is below the water table of the worker gaining water
%                 
%                 [ ind_cT, missingHeightInCell ] = getWTsoilcTindex(soilGrid_workerj,altitude_j,waterpotWindex );% find the index of the cell of worker loosing water matching with the waterTable elevation of worker gaining water
%                 [ ind_wT, ~ ] = getWTsoilcTindex(soilGrid_workerj,altitude_j,PARA.ensemble.waterTable(j) );
%                 if ind_cT==ind_wT;
%                     fprintf('WARNING from getWaterFluxes2Workers : waterTable of workers 1 and 2 in the same gridcell\n')
%                     PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
%                     PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
%                     return
%                 end
%                 cellresidual_water=Wc_workerj(ind_cT)*missingHeightInCell; %err here
%                 availableWater_li(ind_cT:end)=0; % Stop accounting for water after the watertable of the worker gaining water
%                 
%             end
%             
%             availableWater=sum((Wc_workerj(availableWater_li)-PARA.soil.fieldCapacity).*K_deltaSoil_workerj(availableWater_li))+cellresidual_water;
%             
%             lateralFlux_inWaterHeight=DarcyFlux;% min(MaxWaterHeight,availableWater);
%             PARA.ensemble.lateralFlux_inWaterHeight(j,index)=-lateralFlux_inWaterHeight; % lateralFlux_inWaterHeight is a non symetric matrix where mat(worker1,worker2) givesthe water height between both but scale to worker 1
%             PARA.ensemble.lateralFlux_inWaterHeight(index,j)=(PARA.ensemble.weight(j)/PARA.ensemble.weight(index)).*lateralFlux_inWaterHeight;
%             
            water_flux_j=DarcyFlux;
        elseif (waterpotWindex > waterpotWj  && hasWater_index==1) % Current worker is loosing water

            % Calculate maximum of the exchange water volume
            DeltaH= waterpotWindex - waterpotWj;
            contact_height =waterpotWindex - nanmax( [ waterpotWj, inf_index ] );
            Distance=sqrt(DeltaH^2 + PARA.ensemble.hydraulicDistance(j,index)^2);
            %section=(PARA.ensemble.waterTable(index)-GRID.soil.soilGrid(PARA.ensemble.bottomBucketSoilcTIndex(index)+1)) * PARA.ensemble.contactLength(j,index);
            %Jan: see question above
            section=contact_height .* PARA.ensemble.hydraulic_contact_length(j,index);
            DarcyFlux=PARA.ensemble.hydraulic_conductivity(j,index) * (DeltaH/Distance) * section / PARA.ensemble.area(index);
            %MaxWaterHeight=DarcyFlux * PARA.technical.syncTimeStep;
            
            
            
%             % Calculate available water volume that can be lost
%             availableWater_li=wc > PARA.soil.fieldCapacity; % Find where there is water above field capacity
%             cellresidual_water=0;
%             
%             if (altitude_index-GRID.soil.soilGrid(PARA.ensemble.bottomBucketSoilcTIndex(index)+1))> waterpotWj; % Bucket bottom of the worker loosing water is higher than water table of the worker gaining water
%                 
%                 availableWater_li(PARA.ensemble.bottomBucketSoilcTIndex(index)+1:end)=0; % Stop accounting for water after the bottom of the bucket
%                 
%             else % Bucket bottom of the worker loosing water is below the water table of the worker gaining water
%                 
%                 [ ind_cT, missingHeightInCell ] = getWTsoilcTindex(GRID.soil.soilGrid,altitude_index,waterpotWj ); % find the index of the cell of worker loosing water matching with the waterTable elevation of worker gaining water
%                 [ ind_wT, ~ ] = getWTsoilcTindex(GRID.soil.soilGrid,altitude_index,PARA.ensemble.waterTable(index) );
%                 if ind_cT==ind_wT;
%                     fprintf('WARNING from getWaterFluxes2Workers : waterTable of workers 1 and 2 in the same gridcell\n')
%                     PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
%                     PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
%                     return
%                 end
%                 cellresidual_water=missingHeightInCell*wc(ind_cT); % err here
%                 availableWater_li(ind_cT:end)=0; % Stop accounting for water after the watertable of the worker gaining water
%                 
%             end
%             
%             K_deltaSoil=GRID.general.K_delta(GRID.soil.cT_domain);
%             availableWater=sum((wc(availableWater_li)-PARA.soil.fieldCapacity).*K_deltaSoil(availableWater_li))+cellresidual_water;
%             lateralFlux_inWaterHeight=min(MaxWaterHeight,availableWater);
            % Jan: for now: take maximum flux
            
            water_flux_j=-DarcyFlux;
            
%             lateralFlux_inWaterHeight=DarcyFlux;
%             PARA.ensemble.lateralFlux_inWaterHeight(j,index)=(PARA.ensemble.weight(index)/PARA.ensemble.weight(j)).*lateralFlux_inWaterHeight;
%             PARA.ensemble.lateralFlux_inWaterHeight(index,j)=-lateralFlux_inWaterHeight;
        else % same water table, no lateralFlux
            water_flux_j=0;
%             PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
%             PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
        end
        
    else % No possible movement of water, because of snow or suficial freezing
%         PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
%         PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
        water_flux_j=0;
    end
else  
%     % No connection between workers
%     PARA.ensemble.lateralFlux_inWaterHeight(j,index)=NaN;
%     PARA.ensemble.lateralFlux_inWaterHeight(index,j)=NaN;
    water_flux_j = NaN;
end

%water_fluxes = PARA.ensemble.lateralFlux_inWaterHeight(index);


end
