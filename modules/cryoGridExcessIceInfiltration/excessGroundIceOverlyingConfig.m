function [ status ] = excessGroundIceOverlyingConfig( T, wc, GRID )
% Find the start cell based on the mobile water for the excess ground ice
% routine. Modification from Leo, January 19 to cope for the bugs related
% to the melt of excess ice due to lateral contact when the overlying cells
% are frozen

natPor=GRID.general.K_delta(GRID.soil.cT_domain).*GRID.soil.cT_natPor;
water=GRID.general.K_delta(GRID.soil.cT_domain).*wc;
mobileWater = double(T(GRID.soil.cT_domain)>0) .* (water-natPor) .* double(water>natPor); %
[startCell, ~]= LayerIndex(mobileWater~=0);

Tsoil=T(GRID.soil.cT_domain);

if sum( Tsoil(1:startCell) < 0) > 0;
    status='frozen';
else
    status='unfrozen';
end

end

