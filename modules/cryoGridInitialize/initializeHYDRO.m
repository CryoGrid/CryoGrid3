function HYDRO = initializeHYDRO()

% NEW variable names
% storage
HYDRO.dW_soil = 0;
HYDRO.dW_snow = 0;
% precipitation
HYDRO.dp_rain=0;
HYDRO.dp_snow=0; % SWE
% evapotranspiration and sublimation
HYDRO.de=0;
HYDRO.ds=0;
% runoff
HYDRO.dr_surface=0;
HYDRO.dr_subsurface=0;
HYDRO.dr_snowmelt=0;
HYDRO.dr_excessSnow=0;
HYDRO.dr_rain=0;  % this is only rain on frozen ground

end
% OLD variable names:
% HYDRO.rain = 0;
% HYDRO.snow = 0;
% HYDRO.sublimation = 0;
% HYDRO.condensation = 0;
% HYDRO.evapotranspiration = 0;
% HYDRO.surfaceRunoffBucket = 0;
% HYDRO.surfaceRunoffRain = 0;
% HYDRO.surfaceRunoffSnow = 0;
% HYDRO.surfaceRunoffLastCell = 0;
% HYDRO.surfaceRunoffGroundIce = 0;
% HYDRO.lateralFluxBucket = 0;
% HYDRO.excessSnow = 0;
% HYDRO.initialSnow = 0;
