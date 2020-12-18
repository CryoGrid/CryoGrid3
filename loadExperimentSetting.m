function PARA = loadExperimentSetting( PARA )

    PARA.Exp.Case='1Tile';

    PARA.Exp.SetLoc = 'Prudhoe';        % specify location,  current choices:  Prudhoe, Drew, HapVal
    PARA.Exp.Scen='RCP85';                       % specify scenario, current choices: CTR, RCP85
    PARA.Exp.forcing.filename='GFDLbiascor_Prudhoe_RCP85_1975-2095_SinS.mat';

    PARA.modules.xice=0;            % true if thaw subsidence is enabled

    PARA.Exp.zOSL = 0.3;               % depth of organic surface layer [m]
    PARA.snow.rho_snow=250;
    PARA.soil.relative_maxWater=0.;   % relative Max Water (max depth water table above soil) [m]
    PARA.snow.relative_maxSnow=0.40;     % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.Tinitial=[-2 10; 0 5; 0.5 0; 10 -7; 20 -8.6; 600 -1; 1100 11.5]; % temperature profile for Deadhorse
    %Lachenbruch JGR 82, borehole Prudhoe Bay T at 600m -1, below larger gradient, Qgeo=0.05W/m^2 and k=2.0 W/Km ( =(0.3*sqrt(0.57)+0.7*sqrt(3.0))^2, results in ~2.5 K/100m temperature increase

    PARA.IS.RoadOrientation = 0;   % Road Orientation 0: East/West facing -> use standard SW forcing; 1: South facing -> use adapted SW forcing
   
    PARA.Exp.ExpSet = 'Qe_old_road';
    PARA.Exp.OutDir = './Test/';
end


