function waterTable = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table;
% modified such that rounding errors when wc=fieldCapacity are prevented

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
    porosity=GRID.soil.cT_actPor;%(1 - GRID.soil.cT_mineral - GRID.soil.cT_organic);
    i=1;
    eps=1e-9; % this is to prevent rounding errors which lead to wrong water table --> if water>=porosity-eps this is considered saturated
    while water(i)<porosity(i)-eps && T(i)>0
        waterTable = waterTable - K_delta(i);
        i=i+1;
    end
    if i==1
        if T(i)<=0
            waterTable = NaN;
        end
    else
        % check cell above
        if water(i-1)>fieldCapacity(i-1)
            Add = K_delta(i-1) .* ( (water(i-1)-fieldCapacity(i-1)) ./ (porosity(i-1)-fieldCapacity(i-1)) );
            waterTable = waterTable + Add;
        elseif T(i)<=0
            waterTable=NaN;
        end
    end
end