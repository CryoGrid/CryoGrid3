function PARA = loadExperimentSetting( PARA )

    PARA.Exp.Case='1Tile';
    PARA.Exp.SetLoc ='Potsdam';        % specify location,  current choices:  Prudhoe, Drew, HapVal
    PARA.Exp.Scen='ERA5';                       % specify scenario, current choices: CTR, RCP85
    PARA.Exp.forcing.filename='ERA5_Telegrafenberg_1979-2019'
    
    PARA.modules.xice=0;            % true if thaw subsidence is enabled
    PARA.Exp.zOSL = 0.2;               % depth of organic surface layer [m] - not needed in this setting!
    PARA.snow.rho_snow = 250;
    PARA.soil.relative_maxWater = 0.;   % relative Max Water (max depth water table above soil) [m]
   % PARA.snow.relative_maxSnow=0.40;     % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.snow.relative_maxSnow = 2.0;     % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.Tinitial=[-2 25; 0 25; 1 20; 2 16; 4 10; 12 10;  100 13; 1000 40]; %  Potsam timeseries for June & assumption geothermal gradient 3°C/100m 
    PARA.IS.RoadOrientation = 0;   % Road Orientation 0: East/West facing -> use standard SW forcing; 1: South facing -> use adapted SW forcing

    PARA.Exp.ExpSet = 'Test';
    PARA.Exp.OutDir = '/data/permarisk/CryoGrid3/Runs_Potsdam/'
end


