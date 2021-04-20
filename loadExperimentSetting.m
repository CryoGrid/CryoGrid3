function PARA = loadExperimentSetting( PARA )  % vulnerable case  (conservative case: set PARA.modules.xice=0 and PARA.IS.RoadOrientation = 0 )

    PARA.Exp.Case='GravelRoad'
    PARA.Exp.SetLoc='Prudhoe';
    PARA.Exp.Scen='RCP85';   
    PARA.Exp.forcing.filename='GFDLbiascor_Prudhoe_RCP85_1975-2095_SinS.mat'; % forcing file

    PARA.modules.xice=1;            % true if thaw subsidence is enabled

    PARA.Exp.zOSL = 0.3;               % depth of organic surface layer [m]
    PARA.snow.rho_snow=250;
    PARA.soil.relative_maxWater=0.;   % relative Max Water (max depth water table above soil) [m]
    PARA.snow.relative_maxSnow=0.40;     % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.Tinitial=[-2 10; 0 5; 0.5 0; 10 -7; 20 -8.6; 600 -1; 1100 11.5]; % temperature profile for Deadhorse
    %Lachenbruch JGR 82, borehole Prudhoe Bay T at 600m -1, below larger gradient, Qgeo=0.05W/m^2 and k=2.0 W/Km ( =(0.3*sqrt(0.57)+0.7*sqrt(3.0))^2, results in ~2.5 K/100m temperature increase

    PARA.IS.RoadOrientation = 1;   % Road Orientation 0: East/West facing -> use standard SW forcing; 1: South facing -> use adapted SW forcin
   
    PARA.Exp.ExpSet = 'exW2mm_snow4';
    PARA.Exp.Dir = 'Runs_submission/';    
end


