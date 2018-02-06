function c_temp = capacityUnfrozen(wc,GRID,PARA)

c_w = PARA.constants.c_w; % 4.2*10^6; %[J/m�K]
c_o = PARA.constants.c_o; % 2.5*10^6; %[J/m�K]
c_m = PARA.constants.c_m; % 2*10^6; %[J/m�K]

c_temp = (GRID.soil.cT_mineral.*c_m + GRID.soil.cT_organic.*c_o + wc.*c_w);
     
     
% adjust for free water
freeWater_domain = GRID.soil.cT_mineral+GRID.soil.cT_organic<1e-6; % cells without soil matrix material
c_temp(freeWater_domain) = c_w; % assume pure water for cells which consist partly of air and water

% adjust for low soil matrix
lowMinOrg_domain = GRID.soil.cT_mineral+GRID.soil.cT_organic>=1e-6 & ~GRID.soil.excessGroundIce & ( GRID.soil.cT_actPor > GRID.soil.cT_natPor );   % cells with lower soil matrix material than 1-natPor

if sum(lowMinOrg_domain)>0

    water = wc;
    mineral = GRID.soil.cT_mineral;
    organic = GRID.soil.cT_organic;
    matrix = mineral + organic;
    mineral(lowMinOrg_domain) = mineral(lowMinOrg_domain) ./ matrix(lowMinOrg_domain) .* (1 - GRID.soil.cT_natPor(lowMinOrg_domain)) ;
    organic(lowMinOrg_domain) = organic(lowMinOrg_domain) ./ matrix(lowMinOrg_domain) .* (1 - GRID.soil.cT_natPor(lowMinOrg_domain)) ;
    water(lowMinOrg_domain) = min( water(lowMinOrg_domain), GRID.soil.cT_natPor(lowMinOrg_domain) ) ;
    c_temp(lowMinOrg_domain) = water(lowMinOrg_domain) .* c_w + mineral(lowMinOrg_domain) .* c_m + organic(lowMinOrg_domain).* c_o;
end
