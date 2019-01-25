function [GRID, PARA, meltwaterGroundIce, wc] = excessGroundIceThawLateral(T, wc, GRID, PARA)
% Excess ground ice thaw function when the melting is occuring because of
% the lateral fluxes and the overlying cells above the mobile water are
% frozen. Then, the water is removed and the cell shrinked without vertical
% motion of water or soil material. NOTE : for this to work, initial cells
% need to be big enough so that when the excess ice is melted, they are
% still bigger than 2 cm.

% Get the startCell that can shrink
natPor=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_natPor;
water=GRID.general.K_delta(GRID.soil.cT_domain).*wc;
mobileWater = double(T(GRID.soil.cT_domain)>0) .* (water-natPor) .* double(water>natPor); %
[startCell, ~]= LayerIndex(mobileWater~=0);

% Calulate its skrinked height
Kdelta = GRID.general.K_delta(GRID.soil.cT_domain);
finalThick = Kdelta(startCell) * (GRID.soil.cT_mineral(startCell) + GRID.soil.cT_organic(startCell)) / (1-natPor(startCell));
deltaz = finalThick - Kdelta(startCell);
PARA.location.shrinkage = PARA.location.shrinkage + abs(deltaz); 

% Update GRID values
GRID.general.K_grid(GRID.soil.cT_domain_ub + startCell - 1 + 1 : end) = ...
    GRID.general.K_grid(GRID.soil.cT_domain_ub + startCell - 1 + 1 : end) - abs(deltaz);

GRID.general.K_delta(GRID.soil.cT_domain_ub + startCell - 1) = finalThick;

GRID.general.cT_grid(GRID.soil.cT_domain_ub + startCell - 1) = ...
    GRID.general.cT_grid(GRID.soil.cT_domain_ub + startCell - 1) - deltaz/2;

GRID.general.cT_grid(GRID.soil.cT_domain_ub + startCell - 1 + 1 : end) = ...
    GRID.general.cT_grid(GRID.soil.cT_domain_ub + startCell - 1 +1 :end) - deltaz;

GRID.general.cT_delta(GRID.soil.cT_domain_ub + startCell - 2) = ...
    GRID.general.K_delta(GRID.soil.cT_domain_ub + startCell - 2)/2 + finalThick/2; 

GRID.general.cT_delta(GRID.soil.cT_domain_ub + startCell - 1) = ...
    GRID.general.K_delta(GRID.soil.cT_domain_ub + startCell - 1)/2 + finalThick/2;    
% Check calc
cT_delta_check=(-GRID.general.cT_grid(1:end-1,1)+GRID.general.cT_grid(2:end,1));
K_delta_check=(-GRID.general.K_grid(1:end-1,1)+GRID.general.K_grid(2:end,1));
if cT_delta_check ~= GRID.general.cT_delta || K_delta_check ~= GRID.general.K_delta;
    error('excessGroundIceThawLateral : Wrong calculation of the delta grids ') 
end

GRID.soil.cT_mineral(startCell) = GRID.soil.cT_mineral(startCell) * (Kdelta(startCell) / finalThick);
GRID.soil.cT_organic(startCell) = GRID.soil.cT_organic(startCell) * (Kdelta(startCell) / finalThick);
GRID.soil.cT_actPor(startCell)  = 1 - GRID.soil.cT_mineral(startCell) - GRID.soil.cT_organic(startCell);

% Water content
meltwaterGroundIce = water(startCell) - GRID.soil.cT_actPor(startCell);
wc(startCell)      = GRID.soil.cT_actPor(startCell);
warning('xice - reinitializing LUT - modified soil composition')
GRID.soil.cT_water = wc;
GRID = initializeSoilThermalProperties(GRID, PARA);

% Update altitudes
PARA.location.altitude = getAltitude( PARA, GRID );
PARA.location.soil_altitude = getSoilAltitude(PARA, GRID);
PARA.location.surface_altitude = getSurfaceAltitude(PARA, GRID);
[PARA.location.infiltration_altitude, PARA.location.bottomBucketSoilcTIndex] = getInfiltrationAltitude( PARA, GRID, T);
[PARA.location.water_table_altitude, GRID.soil.flag] = getWaterTableAltitudeFC(T, wc, GRID, PARA);
PARA.soil.infiltration_limit_altitude = PARA.location.soil_altitude - PARA.soil.infiltration_limit_depth;

PARA.location.absolute_maxWater_altitude = PARA.location.absolute_maxWater_altitude - abs(deltaz);
PARA.location.absolute_maxSnow_altitude = PARA.location.absolute_maxSnow_altitude - abs(deltaz);


end