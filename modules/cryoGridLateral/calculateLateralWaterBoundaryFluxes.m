<<<<<<< HEAD
function [ boundary_water_fluxes ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T)
% Function that calculates the lateral water fluxes created by the boundary
% conditions of the worker

boundary_water_fluxes=0;
=======
function [ boundary_water_flux ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T)
% Function that calculates the lateral water fluxes created by the boundary
% conditions of the worker

boundary_water_flux=0;
>>>>>>> origin/xice_mpi_polygon_TC

if double( T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub) )==1; % conditions ok for water fluxes
    
    if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1; % worker as the DarcyReservoir boundary condition
        wt=PARA.ensemble.water_table_altitude(labindex);
        inf_altitude = PARA.ensemble.infiltration_altitude(labindex);
        [waterpot, hasWater] = nanmax([wt, inf_altitude] );
        Darcy_elevation=PARA.ensemble.boundaryCondition(labindex).parameters.elevation;
        Darcy_fluxFactor=PARA.ensemble.boundaryCondition(labindex).parameters.fluxFactor;
        DeltaH=abs(waterpot - Darcy_elevation);
        DarcyFlux= Darcy_fluxFactor * DeltaH * DeltaH; % DeltaH is multiplied twice, once as a pressure gradient and the second as the height of the section through which the flux is going
        waterHeight_change=DarcyFlux * PARA.technical.syncTimeStep *24 *3600 / PARA.ensemble.area(labindex); % syncTimeStep in days, Darcy flux in m3/sec
        
        if (waterpot > PARA.ensemble.boundaryCondition(labindex).parameters.elevation && hasWater==1) % worker is loosing water
            
<<<<<<< HEAD
            boundary_water_fluxes=-waterHeight_change;
            
        elseif waterpot < PARA.ensemble.boundaryCondition(labindex).parameters.elevation;
            
            boundary_water_fluxes=waterHeight_change;
            
        else
            
            boundary_water_fluxes=0;
=======
            boundary_water_flux=-waterHeight_change;
            
        elseif waterpot < PARA.ensemble.boundaryCondition(labindex).parameters.elevation;
            
            boundary_water_flux=waterHeight_change;
            
        else
            
            boundary_water_flux=0;
>>>>>>> origin/xice_mpi_polygon_TC
            
        end
        
    end
    
end

end

