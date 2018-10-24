<<<<<<< HEAD
function [waterTable, flag_str] = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table

% Define some variables
T=T(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
fieldC=PARA.soil.fieldCapacity;
flag_str=GRID.soil.flag;
=======
function waterTable = getWaterTableAltitudeFC(T, wc, GRID, PARA)
% Function that find the elevation of the water table;
% modified such that rounding errors when wc=fieldCapacity are prevented

% Define some variables
T=T(GRID.soil.cT_domain);
>>>>>>> origin/xice_mpi_polygon_TC

% Check for favourable conditions for water exchanges and browse through
% the water column
if ~isempty(GRID.snow.cT_domain_ub) || T(1)<0
    waterTable=NaN;
else
<<<<<<< HEAD
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
=======
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

% old implementation of lower part --> led to rounding errors when
% wc=fieldCap
%     while water(i)<=fieldCapacity(i) && T(i)>0
%         waterTable = waterTable - K_delta(i);
%         i=i+1;
%     end
%     
%     if T(i)<=0;
%         waterTable=NaN;
%         % fprintf('Dry to permafrost\n')% Loop was stopped by the temperature condition so no water table
%         
%     elseif water(i)<porosity(i) 
%         Substract=K_delta(i)*( 1 - ( (water(i)-fieldCapacity(i)) / (porosity(i)-fieldCapacity(i)) ));
%         waterTable = waterTable - Substract;
%     end
% end
>>>>>>> origin/xice_mpi_polygon_TC
