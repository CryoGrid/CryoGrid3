function [ boundary_water_flux ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T)
% Function that calculates the lateral water fluxes created by the boundary
% conditions of the worker

boundary_water_flux=0;

if double( T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub) )==1 % conditions ok for water fluxes
    
    if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1 || strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoirNoInflow')==1 % worker as the DarcyReservoir boundary condition
        wt=PARA.ensemble.water_table_altitude(labindex);
        inf_altitude = PARA.ensemble.infiltration_altitude(labindex);
        [waterpot, hasWater] = nanmax([wt, inf_altitude] );
        Darcy_elevation=PARA.ensemble.boundaryCondition(labindex).parameters.elevation;
        Darcy_fluxFactor=PARA.ensemble.boundaryCondition(labindex).parameters.fluxFactor;
        DeltaH=abs(waterpot - Darcy_elevation);
        %contact_height = max( [ waterpot, Darcy_elevation ] ) -  inf_altitude;
        contact_height = abs( waterpot - max( [ Darcy_elevation, inf_altitude ] ) );
        
        DarcyFlux= Darcy_fluxFactor * DeltaH * contact_height; % DeltaH is multiplied twice, once as a pressure gradient and the second as the height of the section through which the flux is going
        waterHeight_change=DarcyFlux * PARA.technical.syncTimeStep *24 *3600 / PARA.ensemble.area(labindex); % syncTimeStep in days, Darcy flux in m3/sec
        
        if (waterpot > Darcy_elevation && hasWater==1) % worker is loosing water
            
            boundary_water_flux=-waterHeight_change;
            
        elseif waterpot < Darcy_elevation % worker is gaining water
            
            if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1
            
                boundary_water_flux=waterHeight_change;
                
            else % noInflow condition
                
                boundary_water_flux=0;
                
            end
            
        else % no gradient
            
            boundary_water_flux=0;
            
        end
        
    elseif strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoirDynamicConductivity')==1
        wt=PARA.ensemble.water_table_altitude(labindex);
        inf_altitude = PARA.ensemble.infiltration_altitude(labindex);
        [waterpot, hasWater] = nanmax([wt, inf_altitude] );
        Darcy_elevation=PARA.ensemble.boundaryCondition(labindex).parameters.elevation;
        % determine conductivity to reservoir based on "amount of degradation"
        Darcy_fluxFactorMax=PARA.ensemble.boundaryCondition(labindex).parameters.fluxFactorMax;
        Darcy_fluxFactorMin=PARA.ensemble.boundaryCondition(labindex).parameters.fluxFactorMax;
        Darcy_elevationMax=PARA.ensemble.boundaryCondition(labindex).parameters.elevationMax;
        Darcy_elevationMin=PARA.ensemble.boundaryCondition(labindex).parameters.elevationMin;
        
        elevation_soil = PARA.ensemble.soil_altitude(labindex);
        
        if elevation_soil > Darcy_elevationMin          % no degradation --> low conductivity to reservoir
            Darcy_fluxFactor = Darcy_fluxFactorMin;
        elseif elevation_soil < Darcy_elevationMax      % high degradation --> high conductivity to reservoir
            Darcy_fluxFactor = Darcy_fluxFactorMax;
        else                                            % intermediate degradation --> linear interpolation
            Darcy_fluxFactor = Darcy_fluxFactorMin + (Darcy_fluxFactorMax - Darcy_fluxFactorMin) .* ( elevation_soil - Darcy_elevationMin ) ./ ( Darcy_elevationMax - Darcy_elevationMin ) ;
        end
       
        
        
        DeltaH=abs(waterpot - Darcy_elevation);
        DarcyFlux= Darcy_fluxFactor * DeltaH * DeltaH; % DeltaH is multiplied twice, once as a pressure gradient and the second as the height of the section through which the flux is going
        waterHeight_change=DarcyFlux * PARA.technical.syncTimeStep *24 *3600 / PARA.ensemble.area(labindex); % syncTimeStep in days, Darcy flux in m3/sec
        
        if (waterpot > PARA.ensemble.boundaryCondition(labindex).parameters.elevation && hasWater==1) % worker is loosing water
            
            boundary_water_flux=-waterHeight_change;
            
        elseif waterpot < PARA.ensemble.boundaryCondition(labindex).parameters.elevation
            
            boundary_water_flux=waterHeight_change;
            
        else
            
            boundary_water_flux=0;
            
        end
        
    end
    
end

end

