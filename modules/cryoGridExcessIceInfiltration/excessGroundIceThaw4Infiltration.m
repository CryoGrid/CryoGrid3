function [GRID, meltwaterGroundIce, wc]=excessGroundIceThaw4Infiltration(T, wc, GRID, PARA)

meltwaterGroundIce=0; % in [m]

%calculates amounts of soil constituents in [m]
mineral=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_mineral;  
organic=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_organic;
natPor=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_natPor;
actPor=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_actPor;

% modification for infiltration
water=GRID.general.K_delta(GRID.soil.cT_domain).*wc;

K_delta=GRID.general.K_delta(GRID.soil.cT_domain);

mobileWater = double(T(GRID.soil.cT_domain)>0) .* (water-natPor) .* double(water>natPor);
[startCell ~]= LayerIndex(mobileWater~=0);

%move solids down
for i=startCell:-1:1
    F_solid_down=K_delta(i)-mineral(i)-organic(i)-natPor(i);
    j=i-1;
    while j>0 && F_solid_down>0
        mineralDown = min(mineral(j), mineral(j)./(mineral(j)+organic(j)).*F_solid_down);
        organicDown = min(organic(j), organic(j)./(mineral(j)+organic(j)).*F_solid_down);
        mineral(i)=mineral(i)+mineralDown;
        organic(i)=organic(i)+organicDown;
        mineral(j)=mineral(j)-mineralDown;
        organic(j)=organic(j)-organicDown;
        F_solid_down=F_solid_down-mineralDown-organicDown;
        j=j-1;
    end
end

%adjust the actual porosity
%natPor(1:startCell)=K_delta(1:startCell)-mineral(1:startCell)-organic(1:startCell);
actPor(1:startCell)=K_delta(1:startCell)-mineral(1:startCell)-organic(1:startCell);

%move water up
mobileWater=0;
for i=startCell:-1:1
    totalWater=water(i)+mobileWater;
    mobileWater=totalWater-actPor(i);
    mobileWater=max(0,mobileWater);
    water(i)=totalWater-mobileWater;
end

%clean up grid cells with non-zero+non-unity water content in domains without soil matrix 
mobileWater=0;
for i=1:startCell
    if mineral(i)+organic(i)==0
        mobileWater=mobileWater+water(i);
        water(i)=0;
    end
end
for i=startCell:-1:1
    if mineral(i)+organic(i)==0
        water_temp=min( [ K_delta(i), mobileWater ] );
        mobileWater=mobileWater-water_temp;
        %water(i)=round(water_temp./K_delta(i)).*K_delta(i);   %this violates the water balance, but ensures that no grid cells with partly water and partly air can exist;
        water(i)=water_temp;
        %water_mismatch = water_temp-water(i);
        %GRID.lake.residualWater = GRID.lake.residualWater + water_mismatch;
        %meltwaterGroundIce=meltwaterGroundIce+water_mismatch; % this corrects the violated water balance: if round=floor then water_mismatch>=0 and this is added to runoff
    end
end

wc=water./K_delta;

GRID.soil.cT_mineral=mineral./K_delta;
GRID.soil.cT_organic=organic./K_delta;
GRID.soil.cT_natPor=natPor./K_delta;
GRID.soil.cT_actPor=actPor./K_delta;
GRID.soil.cT_soilType( (GRID.soil.cT_mineral+GRID.soil.cT_organic)<=1e-6 )=1;  %sets sand freeze curve for all water grid cells

soilGRIDsizeOld = sum( GRID.soil.cT_domain );   % to check if LUT update necessary

% remove air cells and mixed air/water cells above water table and adjust the GRID domains
while GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)+wc(1)<=1e-6 || ...      % upper cell filled with pure air
        (GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)<=1e-6 && ( PARA.location.initial_altitude - GRID.general.K_grid(GRID.soil.cT_domain_ub+1) > PARA.location.absolute_maxWater_altitude) )  

    disp('xice - update GRID - removing grid cell ...')
    if wc(1)<=1e-6
        disp('... upper cell wc=0')
    elseif wc(1)==1
        disp('... upper cell wc=1')
    else
        disp('... upper cell 0<wc<1')
    end
    
    meltwaterGroundIce=meltwaterGroundIce+K_delta(1)*wc(1);

    % adjust air and soil domains and boundaries
    GRID.air.cT_domain(GRID.soil.cT_domain_ub)=1;
    GRID.air.K_domain(GRID.soil.K_domain_ub)=1;
    GRID.air.cT_domain_lb=GRID.air.cT_domain_lb+1;
    GRID.air.K_domain_lb=GRID.air.K_domain_lb+1;
    GRID.soil.cT_domain(GRID.soil.cT_domain_ub)=0;
    GRID.soil.K_domain(GRID.soil.K_domain_ub)=0;
    GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub+1;
    GRID.soil.K_domain_ub=GRID.soil.K_domain_ub+1;
    GRID.soil.soilGrid(1)=[];
    
    %%% modification due to infiltration module
    GRID.soil.cT_water(1)=[];
    wc(1)=[];
    %%%
    
    GRID.soil.cT_organic(1)=[];
    GRID.soil.cT_natPor(1)=[];
    GRID.soil.cT_actPor(1)=[];
    GRID.soil.cT_mineral(1)=[];
    GRID.soil.cT_soilType(1)=[];

    GRID.soil.excessGroundIce(1)=[];

end

% check if the uppermost soil cell contains water above water table
if GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)<=1e-6 && ...
                PARA.location.initial_altitude - GRID.general.K_grid(GRID.soil.cT_domain_ub) > PARA.location.absolute_maxWater_altitude
                
    disp('xice - checking upper cell for excess water');

    actualWater = wc(1)*K_delta(1);
    h = PARA.location.absolute_maxWater_altitude - (PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub+1));
    
    if h<0
        warning('xice - h<0. too much water above water table!')
    end

    if actualWater>h
        wc(1)=h./K_delta(1);
        meltwaterGroundIce = meltwaterGroundIce + actualWater-h;
    end

end    
    
soilGRIDsizeNew = sum (GRID.soil.cT_domain );

% update look up tables since soil water contents changed
cellsChanged = soilGRIDsizeNew - soilGRIDsizeOld;
if cellsChanged > 0
    warning('xice - reinitializing LUT - created soil/water cells in xice module');
    GRID.soil.cT_water = wc;
    GRID = initializeSoilThermalProperties(GRID, PARA);   
elseif cellsChanged < 0
    disp('xice - shortening LUT - removed water cell(s)');
    GRID.soil.cT_water(1) = [];
    GRID.soil.cT_frozen(1) = [];
    GRID.soil.cT_thawed(1) = [];
    GRID.soil.K_frozen(1) = [];
    GRID.soil.K_thawed(1) = [];
    GRID.soil.conductivity(1,:) = [];
    GRID.soil.capacity(1,:) = [];
    GRID.soil.liquidWaterContent(1,:) = [];
end

