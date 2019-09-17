function [PARA, GRID] = surfaceCondition(GRID, PARA, T)   

% set surface parameters (albedo, emissivity, roughnesslength, resistance
% to evaporation) according to the actual surface conditions

GRID.lake.unfrozenWaterSurface=false;

%default soil surface 
PARA.surf.albedo  = PARA.soil.albedo;
PARA.surf.epsilon = PARA.soil.epsilon;
PARA.surf.z0      = PARA.soil.z0;
PARA.surf.rs      = PARA.soil.rs;
    
% check if snow cover exists 
if GRID.snow.cT_domain(GRID.air.cT_domain_lb+1)==1
    PARA.surf.albedo  = PARA.snow.albedo;
    PARA.surf.epsilon = PARA.snow.epsilon;
    PARA.surf.z0      = PARA.snow.z0;
    PARA.surf.rs      = PARA.snow.rs;

% check if water surface exists and whether it is frozen

elseif GRID.soil.cT_domain(GRID.air.cT_domain_lb+1)==1 ...
       && GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-6
       
    % upper soil cell is pure water
    if T(GRID.soil.cT_domain_ub)>0 % unfrozen
        GRID.lake.unfrozenWaterSurface = true;
        PARA.surf.albedo  = PARA.water.albedo;
        PARA.surf.epsilon = PARA.water.epsilon;
        PARA.surf.z0      = PARA.water.z0;
        PARA.surf.rs      = PARA.water.rs;
    else %frozen
        PARA.surf.albedo  = PARA.ice.albedo;
        PARA.surf.epsilon = PARA.ice.epsilon;
        PARA.surf.z0      = PARA.ice.z0;
        PARA.surf.rs      = PARA.ice.rs;
    end
end
   