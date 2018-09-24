function OUT = generateOUT()

    % key state variables
    OUT.cryoGrid3=[];
    OUT.water=[];
    OUT.liquidWater=[]; % for distinction between total and liquid water content
    OUT.timestamp=[];
    OUT.TIMESTEP=[];

    % realated to surface energy balance
    OUT.SEB.Lsta=[];
    OUT.SEB.QE=[];
    OUT.SEB.QH=[];
    OUT.SEB.QNET=[];
    OUT.SEB.QG=[];
    OUT.SEB.Tsurf=[];
    OUT.SEB.albedo_stored=[];
    OUT.SEB.Qsurf=[];

    % related to soil
    OUT.soil.topPosition=[];    % relative to initial altitude
    OUT.soil.lakeFloor=[];      % relative to initial altitude
    OUT.soil.soil=cell(0);

    % related to snow
    OUT.snow.outSnow_i=[];
    OUT.snow.outSnow_a=[];
    OUT.snow.outSnow_w=[];
    OUT.snow.topPosition=[];    % relative to initial altitude
    OUT.snow.botPosition=[];    % relative to initial altitude
    
    % derived characteristics and related to geometry
    OUT.location.area = [];
    OUT.location.altitude=[];
    OUT.location.soil_altitude = [];
    OUT.location.surface_altitude=[];
    OUT.location.infiltration_altitude = [];
    OUT.location.water_table_altitude=[];
    % averaged over output interval
    OUT.location.infiltration_altitude_mean = [];
    OUT.location.water_table_altitude_mean=[];
    

    % lateral fluxes
    OUT.lateral.terrain_index_snow=[];
    OUT.lateral.water_fluxes = [];     % vector containing accumulated lateral water fluxes per output interval in [m] to the current worker from all other workers
    OUT.lateral.snow_flux = [];      % accumulated lateral snow fluxes per output interval in [m SWE] to the current worker
    OUT.lateral.dE_tot = [];      % vector containing depth-integrated lateral heat fluxes per output interval in [J/m^2] to the current worker
    OUT.lateral.dE_cell = [];    % matrix containing cell-wise, accumulated lateral heat fluxes in [J/m^3] to the current worker

    % water balance (WB)
    % all flows are defined as positive when they go into the soil/snow column
    % cumulative values per output interval in [mm]
    % storage
    OUT.WB.dW_soil = [];
    OUT.WB.dW_snow = [];
    % precipitation
    OUT.WB.dp_rain=[];
    OUT.WB.dp_snow=[]; % SWE
    % evapotranspiration and sublimation
    OUT.WB.de=[];
    OUT.WB.ds=[];
    % runoff
    OUT.WB.dr_surface=[];
    OUT.WB.dr_external=[];
    OUT.WB.dr_snowmelt=[];
    OUT.WB.dr_excessSnow=[];
    OUT.WB.dr_lateralSnow=[];   % lateral snow flux to other realizations
    OUT.WB.dr_rain=[];          % this is only rain on frozen ground
    OUT.WB.dr_lateralWater=[];       % lateral water flux to other realizations
    OUT.WB.dr_DarcyReservoir=[]; % lateral water flux created by the boundary condition "DarcyReservoir"
    OUT.WB.dr_lateralExcess=[]; % excess water when applying lateral fluxes
    % mismatches (known)
    OUT.WB.dm_lacking =[];      % mismatch term due to lacking water for evapotranspiration 


    % energy balance (EB)
    % accumulated energy fluxes per output time in [ J / m^2 ]
    OUT.EB.Qg = [];         % ground heat flux (positive into ground)
    OUT.EB.Qe = [];         % latent heat flux (positive into ground)
    OUT.EB.Qh = [];         % sensible heat flux (positive into ground)
    OUT.EB.Qnet = [];
    OUT.EB.Qgeo = [];       % geothermal heat flux
    OUT.EB.Qsurf = [];      % heat flux into uppermost grid cell

    OUT.EB.dE_soil_sens = [];
    OUT.EB.dE_soil_lat = [];
    OUT.EB.dE_soil = [];
    OUT.EB.dE_snow_sens = [];
    OUT.EB.dE_snow_lat = [];
    OUT.EB.dE_snow = [];
    
    OUT.EB.Q_lateral = [];
    
    
    % for DEBUGGING   
    OUT.debugging.dE_dt_SEB = [];
    OUT.debugging.dE_dt_cond = [];
    OUT.debugging.K_grid = [];


end