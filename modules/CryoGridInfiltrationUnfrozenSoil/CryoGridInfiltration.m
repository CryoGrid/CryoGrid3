function [wc, GRID, BALANCE] = CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING, BALANCE, lateral_flux_rate)

if isempty(GRID.snow.cT_domain_ub) && T(GRID.soil.cT_domain_ub)>0   %no snow cover and uppermost grid cell unfrozen
    
    % possible contribution from xice meltwater, snowmelt, rain on frozen ground
    residualWater = GRID.lake.residualWater;
    GRID.lake.residualWater=0;
    
    % external flux
    external_flux_rate = PARA.soil.externalWaterFlux;             % in m/day
    BALANCE.water.dr_external = BALANCE.water.dr_external + external_flux_rate.*timestep.*1000;    %in mm
    
    % lateral flux to/from other workers
    lateral_flux_rate = lateral_flux_rate.*3600.*24; % now in m/day
    BALANCE.water.dr_lateral = BALANCE.water.dr_lateral + lateral_flux_rate.*timestep.*1000;
    
    %%% step 1: infiltrate rain and meltwater and external flux through bucket scheme
    % changes due to evapotranspiration and condensation
    dwc_dt=dwc_dt.*timestep.*24.*3600;   %now in m water per grid cell
    BALANCE.water.de = BALANCE.water.de + sum(dwc_dt)*1000; % in mm accumulated over soil column
    
    % changes due to rainfall
    dwc_dt(1)=dwc_dt(1)+FORCING.i.rainfall./1000.*timestep;
    
    % changes due to residual water
    dwc_dt(1)=dwc_dt(1)+residualWater;
    
    % routing of water
    [wc, surface_runoff, lacking_water] = bucketScheme(T, wc, dwc_dt, GRID, PARA, (external_flux_rate+lateral_flux_rate).*timestep);
    
    % consistency check
    if sum( wc<0 )~=0
        warning( 'CryoGridInfiltration - negative water content occured after bucket scheme' );
        %here one could correct the water balance
        
    end
    
    % remove water above water table in case of ponding, e.g. through rain (independent of xice module)
    if GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)<1e-6 && ...
            PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub)>PARA.location.absolute_maxWater_altitude
        
        
        cellSize = GRID.general.K_delta(GRID.soil.cT_domain_ub);
        actualWater = wc(1)*cellSize;
        h = PARA.location.absolute_maxWater_altitude - (PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub+1));
        if h<0
            warning('h<0. too much water above water table!')
        end
        
        if actualWater>h
            disp('infiltration - removing excess water from upper cell');
            wc(1)=h./cellSize;
            surface_runoff = surface_runoff + actualWater-h;
        end
        
        
    end
    
    %%% step 2: update GRID including reomval of excess water above water table and ponding below water table
    [ wc, GRID, surface_runoff ] = updateGRID_infiltration(wc, GRID, PARA, surface_runoff);
    
    
    % store remaining surface runoff
    BALANCE.water.dr_surface = BALANCE.water.dr_surface - surface_runoff*1000; % in [mm]
    BALANCE.water.dm_lacking = BALANCE.water.dm_lacking + lacking_water*1000;
    
end


% step 3: LUT update
%         JAN:recalculate lookup tables when water content of freezing grid cells
%         has changed (infiltrated cells can freeze --> LUT is updated)
if sum(double(wc~=GRID.soil.cT_water & T(GRID.soil.cT_domain)<=0))>0
    disp('infiltration - reinitializing LUT - freezing of infiltrated cell(s)');
    GRID.soil.cT_water = wc;
    GRID = initializeSoilThermalProperties(GRID, PARA);
end
end