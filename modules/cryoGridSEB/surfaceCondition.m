%tsvd function [PARA, GRID] = surfaceCondition(GRID, PARA, T)    % old implementation
function [PARA, GRID] = surfaceCondition(GRID, PARA, T, t, FORCING, SEB)   

% set surface parameters (albedo, emissivity, roughnesslength, resistance to evaporation) according to the actual surface conditions

% GRID.lake.unfrozenWaterSurface=false; %tsvd not needed anymore

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

% lll comment out this block
% % check if water surface exists and whether it is frozen
% 
% elseif GRID.soil.cT_domain(GRID.air.cT_domain_lb+1)==1 ...
%        && GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-6
%        
%     % upper soil cell is pure water
%     if T(GRID.soil.cT_domain_ub)>0 % unfrozen
%         GRID.lake.unfrozenWaterSurface = true;
%         PARA.surf.albedo  = PARA.water.albedo;
%         PARA.surf.epsilon = PARA.water.epsilon;
%         PARA.surf.z0      = PARA.water.z0;
%         PARA.surf.rs      = PARA.water.rs;
%     else %frozen
%         PARA.surf.albedo  = PARA.ice.albedo;
%         PARA.surf.epsilon = PARA.ice.epsilon;
%         PARA.surf.z0      = PARA.ice.z0;
%         PARA.surf.rs      = PARA.ice.rs;
%     end
% end

%tsvd  check if lake exists
elseif  GRID.lake.water.cT_domain(GRID.air.cT_domain_lb+1)==1 % water surface  
    % GRID.lake.unfrozenWaterSurface = true; %tsvd not needed any more
    %note SolarAzEl.m delivers only an approximation of sun position / t must be in UTC 
    [~, sun_elevation] = SolarAzEl(t,PARA.location.latitude,PARA.location.longitude,PARA.location.altitude);
    PARA.water.albedo = waterAlbedo(sun_elevation, FORCING.i.wind);
    PARA.surf.albedo          = PARA.water.albedo;
    PARA.surf.epsilon         = PARA.water.epsilon; 
    [PARA.surf.z0, z0t, z0q]  = flake_roughnessLength(PARA.water.fetch, FORCING.i.wind, SEB.u_star, 0);
    PARA.surf.z0=real(PARA.surf.z0);
    PARA.surf.rs              = PARA.water.rs;

elseif GRID.lake.ice.z_ice>0 %GRID.lake.ice.cT_domain(GRID.air.cT_domain_lb+1)==1 % ice surface
    PARA.surf.albedo          = PARA.ice.albedo;
    PARA.surf.epsilon         = PARA.ice.epsilon;
    [PARA.surf.z0, z0t, z0q]  = flake_roughnessLength(PARA.water.fetch, FORCING.i.wind, SEB.u_star, 1);
    PARA.surf.z0=real(PARA.surf.z0);
    PARA.surf.rs              = PARA.ice.rs;
end
    
end
