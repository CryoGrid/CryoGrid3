function [waterTableAlt, flag_str] = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table;
% modified such that rounding errors when wc=fieldCapacity are prevented

% Define some variables
T=T(GRID.soil.cT_domain);
flag_str=GRID.soil.flag;

% Check for favourable conditions for water exchanges and browse through
% the water column
if ~isempty(GRID.snow.cT_domain_ub) || T(1)<=0
    waterTableAlt=NaN;
else
    K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
    soilType = GRID.soil.cT_soilType;
    fieldCapacity = zeros(size(soilType));
    [~,i_max]=min(abs((PARA.location.altitude - GRID.soil.soilGrid)-PARA.soil.infiltration_limit_altitude));
    for i=1:size(PARA.soil.soilTypes,1)
        fieldCapacity(soilType==i) = PARA.soil.soilTypes( i, 2 );
    end
    waterTableAlt=PARA.location.altitude;
    water=wc;
    porosity=GRID.soil.cT_actPor;%(1 - GRID.soil.cT_mineral - GRID.soil.cT_organic);
    i=1;
    eps=1e-9; % this is to prevent rounding errors which lead to wrong water table --> if water>=porosity-eps this is considered saturated
    while water(i)<porosity(i)-eps && T(i)>0 && i<=i_max;
        waterTableAlt = waterTableAlt - K_delta(i);
        i=i+1;
    end
    
    % check cell above
    if i>1 && (water(i-1)>fieldCapacity(i-1))
        Add = K_delta(i-1) .* ( (water(i-1)-fieldCapacity(i-1)) ./ (porosity(i-1)-fieldCapacity(i-1)) );
        waterTableAlt = waterTableAlt + Add;
    end
    
    % End conditions
    if T(i)<=0; % Loop stopped because of permafrost was hit but no saturation reached
        
        waterTableAlt=NaN;
        if GRID.soil.flag.dry2permafrost==0;
            fprintf('Dry to permafrost\n')% Loop was stopped by the temperature condition so no water table
            GRID.soil.flag.dry2permafrost=1;
        end
        
    elseif i>i_max; % Max infiltration reached but no saturation
        
        waterTableAlt=NaN;
        if GRID.soil.flag.noWTnoPF==0;
            fprintf('No water table and no permafrost down to infiltration depth\n')% Loop was stopped because no water table was found and no permafrost either
            GRID.soil.flag.noWTnoPF=1;
        end
        
    else
        GRID.soil.flag.dry2permafrost=0;
        GRID.soil.flag.noWTnoPF=0;
    end
    
    flag_str=GRID.soil.flag;
    
end

end