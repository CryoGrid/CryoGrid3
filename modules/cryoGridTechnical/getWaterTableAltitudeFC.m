function [waterTable, flag_str] = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table

% Define some variables
T=T(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
fieldC=PARA.soil.fieldCapacity;
flag_str=GRID.soil.flag;

% Check for favourable conditions for water exchanges and browse through
% the water column
if ~isempty(GRID.snow.cT_domain_ub) || T(1)<0
    waterTable=NaN;
else
    waterTable=PARA.location.altitude;
    water=wc;
    porosity=(1 - GRID.soil.cT_mineral - GRID.soil.cT_organic);
    i=1;
    while water(i)<=fieldC && T(i)>0 && i<=length(water);
        waterTable = waterTable - K_delta(i);
        i=i+1;
    end
    
    % Check ending conditions, adjust wtaertable accordingly and display
    % necessary informations
    if T(i)<=0;
        waterTable=NaN;
        if GRID.soil.flag.dry2permafrost==0;
            fprintf('Dry to permafrost\n')% Loop was stopped by the temperature condition so no water table
            GRID.soil.flag.dry2permafrost=1;
        end
    else
        GRID.soil.flag.dry2permafrost=0;
        
        if i>=length(water);
            waterTable=NaN;
            if GRID.soil.flag.noWTnoPF==0;
                fprintf('No water table and no permafrost\n')% Loop was stopped because no water table was found and no permafrost either
                GRID.soil.flag.noWTnoPF=1;
            end
        else
            GRID.soil.flag.noWTnoPF=0;
            
            if water(i)<porosity(i)
                Substract=K_delta(i)*(1-((water(i)-fieldC)/(porosity(i)-fieldC)));
                waterTable = waterTable - Substract;
            end
        end
    end
    flag_str=GRID.soil.flag;
end