function waterTable = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table

% Define some variables
T=T(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
fieldC=PARA.soil.fieldCapacity;

% Check for favourable conditions for water exchanges and browse through
% the water column
if ~isempty(GRID.snow.cT_domain_ub) || T(1)<0
    waterTable=NaN;
else
    waterTable=PARA.location.altitude;
    water=wc;
    porosity=(1 - GRID.soil.cT_mineral - GRID.soil.cT_organic);
    i=1;
    while water(i)<=fieldC && T(i)>0
        waterTable = waterTable - K_delta(i);
        i=i+1;
%         if i>length(water);
%            water
%            T
%         end
    end
    
    if T(i)<=0;
        waterTable=NaN;
        fprintf('Dry to permafrost\n')% Loop was stopped by the temperature condition so no water table
        
    elseif water(i)<porosity(i) 
        Substract=K_delta(i)*(1-((water(i)-fieldC)/(porosity(i)-fieldC)));
        waterTable = waterTable - Substract;
    end
end