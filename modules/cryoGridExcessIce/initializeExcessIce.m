function [GRID,PARA] = initializeExcessIce(GRID,PARA) 

% set the natural porosity that only water in "mobilewaterDomain" is mobile
if isempty(PARA.soil.mobileWaterDomain) 
    GRID.soil.cT_natPor=GRID.soil.cT_water; % why???
else
    mobileWaterDomain = GRID.general.cT_grid(GRID.soil.cT_domain)>=PARA.soil.mobileWaterDomain(1) & GRID.general.cT_grid(GRID.soil.cT_domain)<=PARA.soil.mobileWaterDomain(2);
    GRID.soil.cT_natPor(~mobileWaterDomain)=GRID.soil.cT_water(~mobileWaterDomain); % why???
end

% lower the water table if air is present above the excess ground ice 
GRID.soil.excessGroundIce = GRID.soil.cT_water>GRID.soil.cT_natPor;
firstCellExcessIce=find(GRID.soil.excessGroundIce(:,1)==1, 1, 'first');

if ~isempty(firstCellExcessIce) && firstCellExcessIce>1
    PARA.soil.waterTable=max(PARA.soil.waterTable,...
                             sum((1 - GRID.soil.cT_water(1:firstCellExcessIce-1) - GRID.soil.cT_mineral(1:firstCellExcessIce-1) - GRID.soil.cT_organic(1:firstCellExcessIce-1))...
                             .*GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+firstCellExcessIce-2)));
end