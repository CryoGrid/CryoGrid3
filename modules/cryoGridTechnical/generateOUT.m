function OUT = generateOUT()

    OUT.snow.outSnow_i=[];
    OUT.snow.outSnow_a=[];
    OUT.snow.outSnow_w=[];

    OUT.cryoGrid3=[];
    OUT.water=[];
    OUT.liquidWater=[]; % for distinction between total and liquid water content
    OUT.timestamp=[];
    OUT.TIMESTEP=[];

    %auxiliary for tracking heat fluxes
    OUT.SEB.dE_dt_SEB = [];
    OUT.SEB.dE_dt_cond = [];

    OUT.SEB.Lsta=[];
    OUT.SEB.QE=[];
    OUT.SEB.QH=[];
    OUT.SEB.QNET=[];
    OUT.SEB.QG=[];
    OUT.SEB.Tsurf=[];
    OUT.SEB.albedo_stored=[];
    OUT.SEB.Qsurf=[];

    OUT.soil.topPosition=[];
    OUT.soil.lakeFloor=[];
    OUT.soil.soil=cell(0);

    OUT.snow.topPosition=[];
    OUT.snow.botPosition=[];

    % for DEBUGGING
    OUT.K_grid = [];

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
    OUT.WB.dr_rain=[];  % this is only rain on frozen ground


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

end