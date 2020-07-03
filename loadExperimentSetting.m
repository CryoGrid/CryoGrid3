function PARA = loadExperimentSetting( PARA )

    PARA.Exp.SetLoc = 'Prudhoe';        % specify location,  current choices:  Prudhoe, Drew, HapVal
    PARA.Exp.Scen='RCP85';                       % specify scenario, current choices: CTR, RCP85

    PARA.modules.xice=1;            % true if thaw subsidence is enabled

    PARA.Exp.zOSL = 0.3;               % depth of organic surface layer [m]
    PARA.snow.rho_snow=250;
    PARA.soil.relative_maxWater=0.;   % relative Max Water (max depth water table above soil) [m]
    PARA.snow.relative_maxSnow=0.40;     % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.Tinitial=[-2 10; 0 5; 0.5 0; 10 -7; 20 -8.6; 600 -1; 1100 11.5]; % temperature profile for Deadhorse
    %Lachenbruch JGR 82, borehole Prudhoe Bay T at 600m -1, below larger gradient, Qgeo=0.05W/m^2 and k=2.0 W/Km ( =(0.3*sqrt(0.57)+0.7*sqrt(3.0))^2, results in ~2.5 K/100m temperature increase

    PARA.IS.RoadOrientation = 0;   % Road Orientation 0: East/West facing -> use standard SW forcing; 1: South facing -> use adapted SW forcing

%     spec='exW0005_snow3_GFDLbc_oT1_Sin'; % field capacity peat increased, Drew Point Forcing (instead of Prudhoe Bay)
%     s1=num2str(PARA.Exp.zOSL,'%.2f'); s2=num2str(PARA.soil.relative_maxWater,'%.2f'); s3=num2str(PARA.snow.relative_maxSnow,'%.2f'); s4=num2str(PARA.snow.rho_snow,'%3i'); 
%     PARA.Exp.ExpSet=['p',s1(3:4),'mw',s2(3:4),'ms',s3(3:4),'rs',s4,spec];

   %spec='exW1mm_snow4_SinS';  
   %s1=num2str(PARA.Exp.zOSL,'%.2f'); s2=num2str(PARA.soil.relative_maxWater,'%.2f'); s3=num2str(PARA.snow.relative_maxSnow,'%.2f'); s4=num2str(PARA.snow.rho_snow,'%3i'); 
   % PARA.Exp.ExpSet=['p',s1(3:4),'mw',s2(3:4),'ms',s3(3:4),'rs',s4,spec];
   
    PARA.Exp.ExpSet = 'exW2mm_snow4';
    PARA.Exp.Dir = 'Runs_ERL_submission/';    
end


