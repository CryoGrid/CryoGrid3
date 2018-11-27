function [TEMPORARY, OUT, BALANCE] = sum_up_output_store(t, T, wc, lwc, timestep, TEMPORARY, BALANCE, PARA, GRID, SEB, OUT, FORCING, saveDir, run_number) 

TEMPORARY.timestep_sum=TEMPORARY.timestep_sum+(timestep*24*3600)*timestep;
TEMPORARY.T_sum=TEMPORARY.T_sum+T.*timestep;
TEMPORARY.Qe_sum=TEMPORARY.Qe_sum+SEB.Qe.*timestep;
TEMPORARY.Qh_sum=TEMPORARY.Qh_sum+SEB.Qh.*timestep;
TEMPORARY.Qnet_sum=TEMPORARY.Qnet_sum+SEB.Qnet.*timestep;
TEMPORARY.Qg_sum=TEMPORARY.Qg_sum+SEB.Qg.*timestep;

TEMPORARY.Qsurf_sum = TEMPORARY.Qsurf_sum + SEB.Qsurf * timestep;
TEMPORARY.dE_dt_SEB_sum = TEMPORARY.dE_dt_SEB_sum + SEB.dE_dt_SEB * timestep;
TEMPORARY.dE_dt_cond_sum = TEMPORARY.dE_dt_cond_sum + SEB.dE_dt_cond * timestep;

TEMPORARY.waterTableElevation_sum = TEMPORARY.waterTableElevation_sum + PARA.location.water_table_altitude * timestep;
TEMPORARY.bottomBucketSoilDepth_sum = TEMPORARY.bottomBucketSoilDepth_sum + PARA.location.infiltration_altitude * timestep;
TEMPORARY.water2pool_sum=TEMPORARY.water2pool_sum + GRID.soil.water2pool * timestep;

%----store in output table --------------------------------------------
if  t==TEMPORARY.outputTime
       
    TEMPORARY.counter = TEMPORARY.counter+1;
   
    %average over timespan:
    TEMPORARY.dt_out=t-TEMPORARY.t_last;
    TEMPORARY.t_last=t; 
    
    TEMPORARY.T_out=TEMPORARY.T_sum./TEMPORARY.dt_out;
	
    TEMPORARY.Qh=TEMPORARY.Qh_sum./TEMPORARY.dt_out;
    TEMPORARY.Qe=TEMPORARY.Qe_sum./TEMPORARY.dt_out;
    TEMPORARY.Qnet=TEMPORARY.Qnet_sum./TEMPORARY.dt_out;
    TEMPORARY.Qg=TEMPORARY.Qg_sum./TEMPORARY.dt_out;
    
    TEMPORARY.Qsurf = TEMPORARY.Qsurf_sum ./ TEMPORARY.dt_out;          
    TEMPORARY.dE_dt_SEB=TEMPORARY.dE_dt_SEB_sum./TEMPORARY.dt_out;
    TEMPORARY.dE_dt_cond=TEMPORARY.dE_dt_cond_sum./TEMPORARY.dt_out;

	TEMPORARY.waterTableElevation=TEMPORARY.waterTableElevation_sum ./ TEMPORARY.dt_out; % This is already stored in OUT.location. But is done here to compare if an average is useful
    TEMPORARY.bottomBucketSoilDepth=TEMPORARY.bottomBucketSoilDepth_sum ./ TEMPORARY.dt_out; % This is already stored in OUT.location. But is done here to compare if an average is useful        

    TEMPORARY.water2pool=TEMPORARY.water2pool_sum ./ TEMPORARY.dt_out;    

    TEMPORARY.timestep_out=TEMPORARY.timestep_sum./TEMPORARY.dt_out;             
    
    %reset sum variables
    TEMPORARY.T_sum(:)=0;
    
    TEMPORARY.Qh_sum=0;
    TEMPORARY.Qe_sum=0;
    TEMPORARY.Qnet_sum=0;
    TEMPORARY.Qg_sum=0;
    TEMPORARY.Qsurf_sum = 0;
    TEMPORARY.dE_dt_SEB_sum=0;
    TEMPORARY.dE_dt_cond_sum=0;

    TEMPORARY.waterTableElevation_sum=0;
    TEMPORARY.bottomBucketSoilDepth_sum=0;
    TEMPORARY.water2pool_sum=0;
    
    TEMPORARY.timestep_sum=0;
    
    
    %------ store new values in OUT struct -----------------------------------
    
    % key state variables
    OUT.cryoGrid3=[OUT.cryoGrid3 [NaN(GRID.air.cT_domain_lb,1); TEMPORARY.T_out(GRID.air.cT_domain_lb+1:end,1)]];
    OUT.water=[OUT.water [ NaN( GRID.soil.cT_domain_ub-1,1) ; wc ] ];
    OUT.liquidWater=[OUT.liquidWater [ NaN( GRID.soil.cT_domain_ub-1,1) ; lwc ] ];
    OUT.TIMESTEP=[OUT.TIMESTEP; TEMPORARY.timestep_out];
    OUT.timestamp=[OUT.timestamp; t]; 
    
    % realated to surface energy balance
    OUT.SEB.Lsta=[OUT.SEB.Lsta; mean(SEB.L_star)];
    OUT.SEB.QE=[OUT.SEB.QE; TEMPORARY.Qe];
    OUT.SEB.QH=[OUT.SEB.QH; TEMPORARY.Qh];
    OUT.SEB.QG=[OUT.SEB.QG; TEMPORARY.Qg];
    OUT.SEB.QNET=[OUT.SEB.QNET; TEMPORARY.Qnet];
    OUT.SEB.Tsurf=[OUT.SEB.Tsurf; TEMPORARY.T_out(GRID.air.cT_domain_lb+1)];
    OUT.SEB.albedo_stored=[OUT.SEB.albedo_stored; PARA.surf.albedo];


    % related to soil
    OUT.soil.topPosition=[OUT.soil.topPosition; -GRID.general.K_grid(GRID.soil.cT_domain_ub)];
    if ~isempty( GRID.lake.cT_domain_ub )
        OUT.soil.lakeFloor = [OUT.soil.lakeFloor; - GRID.general.K_grid(GRID.lake.cT_domain_lb+1) ];
    else
        OUT.soil.lakeFloor = [OUT.soil.lakeFloor; NaN];
    end
    OUT.soil.soil{1, size(OUT.soil.soil,2)+1}=[wc GRID.soil.cT_mineral GRID.soil.cT_organic];


    % related to snow
    OUT.snow.outSnow_i=[OUT.snow.outSnow_i GRID.snow.Snow_i];
    OUT.snow.outSnow_a=[OUT.snow.outSnow_a GRID.snow.Snow_a];
    OUT.snow.outSnow_w=[OUT.snow.outSnow_w GRID.snow.Snow_w];
    if ~isempty(GRID.snow.cT_domain_ub)
        TEMPORARY.topPosition=-GRID.general.K_grid(GRID.snow.cT_domain_ub);
        TEMPORARY.botPosition=-GRID.general.K_grid(GRID.snow.cT_domain_lb+1);
    else
        TEMPORARY.topPosition=NaN;
        TEMPORARY.botPosition=NaN;
    end
    OUT.snow.topPosition=[OUT.snow.topPosition; TEMPORARY.topPosition];
    OUT.snow.botPosition=[OUT.snow.botPosition; TEMPORARY.botPosition];


    % derived characteristics and related to geometry
    OUT.location.area = [OUT.location.area; PARA.location.area];
    OUT.location.altitude=[OUT.location.altitude; PARA.location.altitude];
    OUT.location.soil_altitude= [OUT.location.soil_altitude; PARA.location.soil_altitude];
    OUT.location.surface_altitude=[OUT.location.surface_altitude; PARA.location.surface_altitude];
    OUT.location.infiltration_altitude = [OUT.location.infiltration_altitude; PARA.location.infiltration_altitude];
    OUT.location.water_table_altitude=[OUT.location.water_table_altitude; PARA.location.water_table_altitude];
    
    OUT.location.infiltration_altitude_mean = [OUT.location.infiltration_altitude_mean ; TEMPORARY.bottomBucketSoilDepth];
    OUT.location.water_table_altitude_mean = [OUT.location.water_table_altitude_mean ; TEMPORARY.waterTableElevation];
    
    % lateral fluxes
    if PARA.modules.lateral
        %OUT.lateral.terrain_index_snow=[ OUT.lateral.terrain_index_snow; PARA.ensemble.terrain_index_snow ];
        OUT.lateral.snow_scaling = [ OUT.lateral.snow_scaling; PARA.ensemble.snow_scaling ];
        OUT.lateral.water_fluxes = cat(3,OUT.lateral.water_fluxes, BALANCE.water.dr_water_fluxes_out );     % vector containing water fluxes in [m/s] to the current worker
        %OUT.lateral.snow_fluxes = [ OUT.lateral.snow_fluxes; snow_fluxes ];                      % vector containing snow fluxes in [m SWE / s] to the current worker
        %OUT.lateral.heat_fluxes = [ OUT.lateral.heat_fluxes; heat_fluxes ];                      % vector containing depth-integrated heat fluxes in [W/m^2] to the current worker
        %OUT.lateral.snow_flux = [ OUT.lateral.snow_flux;  TEMPORARY.snow_flux_lateral ];      % accumulated lateral snow fluxes per output interval in [m SWE] to the current worker
        OUT.lateral.dE_tot = [ OUT.lateral.dE_tot ; TEMPORARY.dE_tot_lateral];      % vector containing depth-integrated lateral heat fluxes per output interval in [J/m^2] to the current worker
        OUT.lateral.dE_cell = cat( 3, OUT.lateral.dE_cell, TEMPORARY.dE_cell_lateral );    % matrix containing cell-wise, accumulated lateral heat fluxes in [J/m^3] to the current worker
        
        TEMPORARY.snow_flux_lateral = 0 ;
        TEMPORARY.dE_cell_lateral = zeros( length(GRID.general.cT_grid), numlabs );
        TEMPORARY.dE_tot_lateral = zeros( 1, numlabs ) ;
    
    end
    
    % related to carbon
    OUT.carbon.unfrozen_organic_volume_time = [ OUT.carbon.unfrozen_organic_volume_time;  TEMPORARY.unfrozen_organic_volume_time ];
    TEMPORARY.unfrozen_organic_volume_time=0;
    OUT.carbon.unfrozen_organic_volume_time_aerobic = [ OUT.carbon.unfrozen_organic_volume_time_aerobic;  TEMPORARY.unfrozen_organic_volume_time_aerobic ];
    TEMPORARY.unfrozen_organic_volume_time_aerobic=0;    
    OUT.carbon.unfrozen_organic_volume_time_anaerobic = [ OUT.carbon.unfrozen_organic_volume_time_anaerobic;  TEMPORARY.unfrozen_organic_volume_time_anaerobic ];
    TEMPORARY.unfrozen_organic_volume_time_anaerobic=0;    
    % related to excess ice
    OUT.xice.excessIceThawed = [ OUT.xice.excessIceThawed; TEMPORARY.excessIceThawed ] ;
    TEMPORARY.excessIceThawed = 0;
    
    % water balance (WB)
    % all flows are defined as positive when they go into the soil/snow column
    % cumulative values per output interval in [mm]
    % storage
    OUT.WB.dW_soil = [ OUT.WB.dW_soil; BALANCE.water.dW_soil ];
    OUT.WB.dW_snow = [ OUT.WB.dW_snow; BALANCE.water.dW_snow ];
    OUT.WB.water2pool= [OUT.WB.water2pool; TEMPORARY.water2pool * 1000]; % Non additive variable. Times 1000 to have it in mm
    % precipitation
    OUT.WB.dp_rain = [ OUT.WB.dp_rain; BALANCE.water.dp_rain ];
    OUT.WB.dp_snow = [ OUT.WB.dp_snow; BALANCE.water.dp_snow ]; % SWE
    % evapotranspiration and sublimation
    OUT.WB.de = [ OUT.WB.de; BALANCE.water.de ];
    OUT.WB.ds = [ OUT.WB.ds; BALANCE.water.ds ];
    % runoff
    OUT.WB.dr_surface= [ OUT.WB.dr_surface; BALANCE.water.dr_surface ];
    OUT.WB.dr_external = [ OUT.WB.dr_external; BALANCE.water.dr_external ];
    OUT.WB.dr_snowmelt = [ OUT.WB.dr_snowmelt; BALANCE.water.dr_snowmelt ];
    OUT.WB.dr_excessSnow=[ OUT.WB.dr_excessSnow; BALANCE.water.dr_excessSnow ];
 	OUT.WB.dr_lateralSnow=[ OUT.WB.dr_lateralSnow; BALANCE.water.dr_lateralSnow ];  % lateral snow flux to other realizations
    OUT.WB.dr_rain = [ OUT.WB.dr_rain; BALANCE.water.dr_rain ];  					% this is only rain on frozen ground
    OUT.WB.dr_lateralWater = [OUT.WB.dr_lateralWater; BALANCE.water.dr_lateralWater ];			% lateral water flux to other realizations
    OUT.WB.dr_DarcyReservoir = [OUT.WB.dr_DarcyReservoir; BALANCE.water.dr_DarcyReservoir ]; % fluxes related to the Darcy Reservoir boundary condition
    OUT.WB.dr_lateralExcess = [OUT.WB.dr_lateralExcess ; BALANCE.water.dr_lateralExcess ]; % % excess water when applying lateral fluxes
    % mismatch
    OUT.WB.dm_lacking = [OUT.WB.dm_lacking; BALANCE.water.dm_lacking ];
    % set accumulated fluxes in BALANCE.water struct to zero
    % storage
    BALANCE.water.dW_soil = 0;
    BALANCE.water.dW_snow = 0;
    % precipitation
    BALANCE.water.dp_rain=0;
    BALANCE.water.dp_snow=0; % SWE
    % evapotranspiration and sublimation
    BALANCE.water.de=0;
    BALANCE.water.ds=0;
    % runoff
    BALANCE.water.dr_surface=0;
    BALANCE.water.dr_external=0;
    BALANCE.water.dr_snowmelt=0;
    BALANCE.water.dr_excessSnow=0;
	BALANCE.water.dr_lateralSnow=0;
    BALANCE.water.dr_rain=0;
    BALANCE.water.dr_lateralWater=0;
    BALANCE.water.dr_DarcyReservoir=0;
    BALANCE.water.dr_lateralExcess=0;
    BALANCE.water.dr_water_fluxes_out=zeros(numlabs,numlabs);
    
    %mismatch
    BALANCE.water.dm_lacking=0;

    % complete energy balance (EB)
    OUT.EB.Qg = [OUT.EB.Qg; TEMPORARY.Qg];         % ground heat flux (positive into ground)
    OUT.EB.Qe = [OUT.EB.Qe; TEMPORARY.Qe];         % latent heat flux (positive into ground)
    OUT.EB.Qh = [OUT.EB.Qh; TEMPORARY.Qh];         % sensible heat flux (positive into ground)
    OUT.EB.Qnet = [OUT.EB.Qnet; TEMPORARY.Qnet];
    OUT.EB.Qgeo = [OUT.EB.Qgeo; PARA.soil.Qgeo];       % geothermal heat flux
	OUT.EB.dE_soil_sens = [OUT.EB.dE_soil_sens; BALANCE.energy.dE_soil_sens ];
    BALANCE.energy.dE_soil_sens = 0;
    OUT.EB.dE_soil_lat = [OUT.EB.dE_soil_lat; BALANCE.energy.dE_soil_lat ];
    BALANCE.energy.dE_soil_lat = 0;
    OUT.EB.dE_soil = [OUT.EB.dE_soil; BALANCE.energy.dE_soil ];
    BALANCE.energy.dE_soil = 0;
    OUT.EB.dE_snow_sens = [ OUT.EB.dE_snow_sens; BALANCE.energy.dE_snow_sens ];
    BALANCE.energy.dE_snow_sens = 0;
    OUT.EB.dE_snow_lat = [ OUT.EB.dE_snow_lat; BALANCE.energy.dE_snow_lat ];
    BALANCE.energy.dE_snow_lat = 0;
    OUT.EB.dE_snow = [ OUT.EB.dE_snow; BALANCE.energy.dE_snow ];
    BALANCE.energy.dE_snow = 0;
    OUT.EB.Q_lateral = [OUT.EB.Q_lateral, [ BALANCE.energy.Q_lateral ] ];
    BALANCE.energy.Q_lateral = zeros( length(GRID.general.cT_grid) , 1 );
        
    % for DEBUGGING   
    OUT.debugging.dE_dt_SEB = [OUT.debugging.dE_dt_SEB [ TEMPORARY.dE_dt_SEB ] ];
    OUT.debugging.dE_dt_cond = [OUT.debugging.dE_dt_cond [ TEMPORARY.dE_dt_cond ] ];
    OUT.debugging.K_grid = [OUT.debugging.K_grid, GRID.general.K_grid ];


    %------------------------------------------------------------------     
    disp([datestr(now,'yyyy-mm-dd HH:MM:SS'),':  at ',datestr(t), ',  Average timestep: ',  num2str(TEMPORARY.timestep_out), ' seconds'])
  
    TEMPORARY.outputTime=round((TEMPORARY.outputTime+PARA.technical.outputTimestep)./PARA.technical.outputTimestep).*PARA.technical.outputTimestep;
 
    %write output files      
    if  round((t-TEMPORARY.saveTime).*48)==0
        DIAG = diagnose_output_yearly( OUT, PARA, GRID, FORCING );
        iSaveDIAG( [ saveDir '/' run_number '/' run_number '_realization' num2str(labindex) '_diagnostics' datestr(t,'yyyy')  '.mat' ], DIAG);      
        iSaveOUT( [ saveDir '/' run_number '/' run_number '_realization' num2str(labindex) '_output' datestr(t,'yyyy')  '.mat' ], OUT);
        iSaveState( [ saveDir '/' run_number '/' run_number '_realization' num2str(labindex) '_finalState'  datestr(t,'yyyy') '.mat' ], T, wc, t, SEB, PARA, GRID);
        iPlotAltitudes( [ saveDir '/' run_number '/' run_number '_realization' num2str(labindex) '_altitudes' datestr(t,'yyyy') '.png'], OUT, PARA );
        OUT = generateOUT();  
        TEMPORARY.saveTime=datenum(str2num(datestr(t,'yyyy'))+1, str2num(datestr(t,'mm')), str2num(datestr(t,'dd')), str2num(datestr(t,'HH')), str2num(datestr(t,'MM')), 0);
    end
end
end
