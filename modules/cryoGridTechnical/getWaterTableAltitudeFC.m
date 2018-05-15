function waterTable = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table

% Define some variables
T=T(GRID.soil.cT_domain);

% Check for favourable conditions for water exchanges and browse through
% the water column
if ~isempty(GRID.snow.cT_domain_ub) || T(1)<0
    waterTable=NaN;
else
    K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
    soilType = GRID.soil.cT_soilType;
    fieldCapacity = zeros(size(soilType));
    for i=1:size(PARA.soil.soilTypes,1)
        fieldCapacity(soilType==i) = PARA.soil.soilTypes( i, 2 );
    end
    
    waterTable=PARA.location.altitude;
    water=wc;
    porosity=(1 - GRID.soil.cT_mineral - GRID.soil.cT_organic);
    i=1;
    while water(i)<=fieldCapacity(i) && T(i)>0
        waterTable = waterTable - K_delta(i);
        i=i+1;
    end
    
    if T(i)<=0;
        waterTable=NaN;
        % fprintf('Dry to permafrost\n')% Loop was stopped by the temperature condition so no water table
        
    elseif water(i)<porosity(i) 
        Substract=K_delta(i)*( 1 - ( (water(i)-fieldCapacity(i)) / (porosity(i)-fieldCapacity(i)) ));
        waterTable = waterTable - Substract;
    end
end