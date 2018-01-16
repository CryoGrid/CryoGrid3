function [GRID, meltwaterGroundIce, PARA]=excessGroundIceThaw4(T, GRID, PARA)

%disp('rearranging grid cells due to ground ice thaw')

waterLevel=PARA.soil.waterTable;     %remove supersaturation only when there is no snow on top of soil!!!!!!!


mineral=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_mineral;  %calculates amounts in [m]
organic=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_organic;
water=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_water;
natPor=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_natPor;

cT_grid=GRID.general.cT_grid(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);


mobileWater = double(T(GRID.soil.cT_domain)>0) .* (water-natPor) .* double(water>natPor);
[startCell ~]= LayerIndex(mobileWater~=0); %this is faster

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

%adjust the natural porosity
natPor(1:startCell)=K_delta(1:startCell)-mineral(1:startCell)-organic(1:startCell);

%move water up
mobileWater=0;
for i=startCell:-1:1
    totalWater=water(i)+mobileWater;
    mobileWater=totalWater-natPor(i);
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
        water(i)=min(K_delta(i), mobileWater);
        mobileWater=mobileWater-water(i);
        water(i)=round(water(i)./K_delta(i)).*K_delta(i);  %this violates the water balance, but ensures that no grid cells with partly water and partly air can exist;
    end
end


GRID.soil.cT_mineral=mineral./K_delta;
GRID.soil.cT_organic=organic./K_delta;
GRID.soil.cT_water=water./K_delta;
GRID.soil.cT_natPor=natPor./K_delta;
test=mineral+water+organic;
GRID.soil.cT_soilType(water(:,1)==1,1)=1;  %sets sand freeze curve for all water grid cells (comment: find was replaced!!!)

GRID.soil.K_mineral(1)=GRID.soil.cT_mineral(1);
GRID.soil.K_mineral(2:startCell+1)=(GRID.soil.cT_mineral(2:startCell+1)+GRID.soil.cT_mineral(1:startCell))/2 ;
GRID.soil.K_organic(1)=GRID.soil.cT_organic(1);
GRID.soil.K_organic(2:startCell+1)=(GRID.soil.cT_organic(2:startCell+1)+GRID.soil.cT_organic(1:startCell))/2 ;
GRID.soil.K_water(1)=GRID.soil.cT_water(1);
GRID.soil.K_water(2:startCell+1)=(GRID.soil.cT_water(2:startCell+1)+GRID.soil.cT_water(1:startCell))/2 ;
% GRID.soil.K_natPor(1)=GRID.soil.cT_natPor(1);
% GRID.soil.K_natPor(2:startCell+1)=(GRID.soil.cT_natPor(2:startCell+1)+GRID.soil.cT_natPor(1:startCell))/2 ;
GRID.soil.K_soilType(1)=GRID.soil.cT_soilType(1);
GRID.soil.K_soilType(2:startCell+1)=round((GRID.soil.cT_soilType(2:startCell+1)+GRID.soil.cT_soilType(1:startCell))/2) ;

meltwaterGroundIce=0;

while (GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)+GRID.soil.cT_water(1)==0) || (GRID.soil.cT_water(1)==1 && GRID.general.K_grid(GRID.soil.K_domain_ub)<waterLevel)   %remove grid cells until the water level is reached (comment: why GRID.general.K_grid(GRID.soil.K_domain_ub)<waterLevel) 
    
    GRID.air.cT_domain(GRID.soil.cT_domain_ub)=1;
    GRID.air.K_domain(GRID.soil.K_domain_ub)=1;
    GRID.air.cT_domain_lb=GRID.air.cT_domain_lb+1;
    GRID.air.K_domain_lb=GRID.air.K_domain_lb+1;
   
    GRID.soil.cT_domain(GRID.soil.cT_domain_ub)=0;
    GRID.soil.K_domain(GRID.soil.K_domain_ub)=0;
    GRID.soil.soilGrid(1)=[];
    %GRID.soil.convectiveDomain(1)=[];
    GRID.soil.cT_water(1)=[];
    GRID.soil.cT_organic(1)=[];
    GRID.soil.cT_natPor(1)=[];
    GRID.soil.cT_mineral(1)=[];
    GRID.soil.cT_soilType(1)=[];
    GRID.soil.K_water(1)=[];
    GRID.soil.K_organic(1)=[];
    GRID.soil.K_mineral(1)=[];
    GRID.soil.K_soilType(1)=[];
    GRID.soil.excessGroundIce(1)=[];
    meltwaterGroundIce=meltwaterGroundIce+GRID.general.K_delta(GRID.soil.cT_domain_ub);
    GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub+1;
    GRID.soil.K_domain_ub=GRID.soil.K_domain_ub+1;
end

GRID = initializeSoilThermalProperties( GRID, PARA );

%reduce water content above the perched water table
%it might be also good to use an external function for this since this is
%only an option in the model
%i=0;
% while     i./startCell < 1-PARA.soil.perchedWaterTable
%     GRID.soil.cT_water(i+1,1) = PARA.soil.saturationAbovePerchedWaterTable.*(1-GRID.soil.cT_mineral(i+1,1)-GRID.soil.cT_organic(i+1,1));
%    
%     GRID.soil.K_water(i+1,1) = PARA.soil.saturationAbovePerchedWaterTable.*(1-GRID.soil.K_mineral(i+1,1)-GRID.soil.K_organic(i+1,1));
%     i=i+1;
% end


% [GRID.soil.cT_frozen,...
%  GRID.soil.cT_thawed,...
%  GRID.soil.K_frozen,...
%  GRID.soil.K_thawed,...
%  GRID.soil.conductivity,...
%  GRID.soil.capacity] = initialize(GRID.soil.cT_water,...
%                                   GRID.soil.cT_mineral,...
%                                   GRID.soil.cT_organic,...
%                                   GRID.soil.cT_soilType,...
%                                   GRID.soil.K_water,...
%                                   GRID.soil.K_mineral,...
%                                   GRID.soil.K_organic,...
%                                   GRID.soil.K_soilType,...
%                                   PARA.technical.arraySizeT,...
%                                   GRID.general.cT_grid(GRID.soil.cT_domain),...
%                                   PARA.soil.kh_bedrock);
