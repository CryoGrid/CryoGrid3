function [ PARA ] = calculateLateralWaterAvailable( PARA,GRID,T, packageWorkerj,index,j, wc )
% Function that calculate the water fluxes between the 2 workers.
% PARA.ensemble should contain the elevations, the water table values, the
% bucketbottom the weight and conductivities of each worcker.
% GRID will contain the soil water content.
% T is the temperature vector of the current worker.
% TWcworkerj has to be a bundle the T profile and the water content of the
% other worker IN THE SOIL, so use he logical indexing before sending

% Check if the 2 workers are connected
if isnan(PARA.ensemble.hydraulic_conductivity(j,index))~=1;
    
    % Prepare data from the other worker
    T_soilub_workerj=packageWorkerj.T_soilub_workerj;
    Wc_workerj=packageWorkerj.Wc_workerj;
    soilGrid_workerj=packageWorkerj.soilGrid_workerj;
    K_deltaSoil_workerj=packageWorkerj.K_deltaSoil_workerj;
    SnowCoverBoolean_workerj=packageWorkerj.SnowCoverBoolean_workerj;
    altitude_index=PARA.ensemble.altitude(index);
    altitude_j=PARA.ensemble.altitude(j);
    
    
    if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0 && T_soilub_workerj>0 && SnowCoverBoolean_workerj==0; % conditions for water exchanges
        
        % Decipher between cases
        [waterpotWindex, i_index] = nanmax([PARA.ensemble.waterTable(index) (PARA.ensemble.altitude(index) - GRID.soil.soilGrid(PARA.ensemble.bottomBucketSoilcTIndex(index)+1))]);
        [waterpotWj, i_j]         = nanmax([PARA.ensemble.waterTable(j) (PARA.ensemble.altitude(j) - soilGrid_workerj(PARA.ensemble.bottomBucketSoilcTIndex(j)+1))]);
        
        if (waterpotWj > waterpotWindex && i_j==1);% Current worker is gaining water
            % fprintf('Worker %1.0f\t has the lower WT\n',index)
            % Calculate the maximum exchanged water volume
            DeltaH=PARA.ensemble.waterTable(j)-waterpotWindex;
            Distance=sqrt(DeltaH^2 + PARA.ensemble.distanceBetweenPoints(j,index)^2);
            section=(PARA.ensemble.waterTable(j)-soilGrid_workerj(PARA.ensemble.bottomBucketSoilcTIndex(j)+1)) * PARA.ensemble.contactLength(j,index);
            DarcyFlux=PARA.ensemble.hydraulic_conductivity(j,index) * (DeltaH/Distance) * section;
            MaxWaterHeight=DarcyFlux * PARA.technical.syncTimeStep / PARA.ensemble.weight(j);
            
            % Calculate available water volume that can be lost
            availableWater_li=Wc_workerj > PARA.soil.fieldCapacity; % Find where there is water above field capacity
            cellresidual_water=0;
            
            if (altitude_j-soilGrid_workerj(PARA.ensemble.bottomBucketSoilcTIndex(j)+1))> waterpotWindex; % Bucket bottom of the worker loosing water is higher than water table of the worker gaining water
                
                availableWater_li(PARA.ensemble.bottomBucketSoilcTIndex(j)+1:end)=0; % Stop accounting for water after the bottom of the bucket
                
            else % Bucket bottom of the worker loosing water is below the water table of the worker gaining water
                
                [ ind_cT, missingHeightInCell ] = getWTsoilcTindex(soilGrid_workerj,altitude_j,waterpotWindex );% find the index of the cell of worker loosing water matching with the waterTable elevation of worker gaining water
                [ ind_wT, ~ ] = getWTsoilcTindex(soilGrid_workerj,altitude_j,PARA.ensemble.waterTable(j) );
                if ind_cT==ind_wT;
                    fprintf('WARNING from getWaterFluxes2Workers : waterTable of workers 1 and 2 in the same gridcell\n')
                    PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
                    PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
                    return
                end
                cellresidual_water=Wc_workerj(ind_cT)*missingHeightInCell; %err here
                availableWater_li(ind_cT:end)=0; % Stop accounting for water after the watertable of the worker gaining water
                
            end
            
            availableWater=sum((Wc_workerj(availableWater_li)-PARA.soil.fieldCapacity).*K_deltaSoil_workerj(availableWater_li))+cellresidual_water;
            lateralFlux_inWaterHeight=min(MaxWaterHeight,availableWater);
            PARA.ensemble.lateralFlux_inWaterHeight(j,index)=-lateralFlux_inWaterHeight; % lateralFlux_inWaterHeight is a non symetric matrix where mat(worker1,worker2) givesthe water height between both but scale to worker 1
            PARA.ensemble.lateralFlux_inWaterHeight(index,j)=(PARA.ensemble.weight(j)/PARA.ensemble.weight(index)).*lateralFlux_inWaterHeight;
            
        elseif (waterpotWindex > waterpotWj  && i_index==1); % Current worker is loosing water
            % fprintf('Worker %1.0f\t has the higher WT\n',index)
            % Calculate maximum of the exchange water volume
            DeltaH=PARA.ensemble.waterTable(index)-waterpotWj;
            Distance=sqrt(DeltaH^2 + PARA.ensemble.distanceBetweenPoints(j,index)^2);
            section=(PARA.ensemble.waterTable(index)-GRID.soil.soilGrid(PARA.ensemble.bottomBucketSoilcTIndex(index)+1)) * PARA.ensemble.contactLength(j,index);
            DarcyFlux=PARA.ensemble.hydraulic_conductivity(j,index) * (DeltaH/Distance) * section;
            MaxWaterHeight=DarcyFlux * PARA.technical.syncTimeStep / PARA.ensemble.weight(index);
            %             if index==2;
            %                 fprintf('Max flux : %3.2e\n',MaxWaterHeight)
            %                 fprintf('Delta H : %3.2f\n', DeltaH)
            %                 [ ind_cT, ~] =getWTsoilcTindex(GRID.general.K_delta(GRID.soil.cT_domain) , PARA.ensemble.altitude(index) , PARA.ensemble.waterTable(index));
            %                 fprintf('index top : %1.0f\n', ind_cT)
            %                 fprintf('index bottom : %1.0f\n',PARA.ensemble.bottomBucketSoilcTIndex(index))
            %             end
            
            % Calculate available water volume that can be lost
            availableWater_li=wc > PARA.soil.fieldCapacity; % Find where there is water above field capacity
            cellresidual_water=0;
            
            if (altitude_index-GRID.soil.soilGrid(PARA.ensemble.bottomBucketSoilcTIndex(index)+1))> waterpotWj; % Bucket bottom of the worker loosing water is higher than water table of the worker gaining water
                
                availableWater_li(PARA.ensemble.bottomBucketSoilcTIndex(index)+1:end)=0; % Stop accounting for water after the bottom of the bucket
                
            else % Bucket bottom of the worker loosing water is below the water table of the worker gaining water
                
                [ ind_cT, missingHeightInCell ] = getWTsoilcTindex(GRID.soil.soilGrid,altitude_index,waterpotWj ); % find the index of the cell of worker loosing water matching with the waterTable elevation of worker gaining water
                [ ind_wT, ~ ] = getWTsoilcTindex(GRID.soil.soilGrid,altitude_index,PARA.ensemble.waterTable(index) );
                if ind_cT==ind_wT;
                    fprintf('WARNING from getWaterFluxes2Workers : waterTable of workers 1 and 2 in the same gridcell\n')
                    PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
                    PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
                    return
                end
                cellresidual_water=missingHeightInCell*wc(ind_cT); % err here
                availableWater_li(ind_cT:end)=0; % Stop accounting for water after the watertable of the worker gaining water
                
            end
            
            K_deltaSoil=GRID.general.K_delta(GRID.soil.cT_domain);
            availableWater=sum((wc(availableWater_li)-PARA.soil.fieldCapacity).*K_deltaSoil(availableWater_li))+cellresidual_water;
            lateralFlux_inWaterHeight=min(MaxWaterHeight,availableWater);
%             if index==2;
%                 fprintf('available flux : %3.2e m\n',availableWater)
%             end
            PARA.ensemble.lateralFlux_inWaterHeight(j,index)=(PARA.ensemble.weight(index)/PARA.ensemble.weight(j)).*lateralFlux_inWaterHeight;
            PARA.ensemble.lateralFlux_inWaterHeight(index,j)=-lateralFlux_inWaterHeight;
        else % same water table, no lateralFlux
            PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
            PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
        end
        
    else % No possible movement of water, because of snow of suficial freezing
        PARA.ensemble.lateralFlux_inWaterHeight(j,index)=0;
        PARA.ensemble.lateralFlux_inWaterHeight(index,j)=0;
    end
    
else % No connection between workers
    PARA.ensemble.lateralFlux_inWaterHeight(j,index)=NaN;
    PARA.ensemble.lateralFlux_inWaterHeight(index,j)=NaN;
end

end

