function k_temp = conductivityUnfrozen(wc,GRID, PARA)

ka = PARA.constants.k_a; %0.025;       %air [Hillel(1982)]
kw = PARA.constants.k_w; %0.57;        %water [Hillel(1982)]
ko = PARA.constants.k_o; %0.25;        %organic [Hillel(1982)]
km = PARA.constants.k_m; %soil.kh_bedrock;     %mineral 

k_freeWater = 5.0; % to ensure fast heat transfer

air=1-wc-GRID.soil.cT_mineral-GRID.soil.cT_organic;

k_temp = (wc.* kw.^0.5 + GRID.soil.cT_mineral.* km.^0.5 + GRID.soil.cT_organic.* ko.^0.5 + air.* ka.^0.5).^2 ;

% adjust for free water
freeWater_domain = GRID.soil.cT_mineral+GRID.soil.cT_organic<1e-6; % cells without soil matrix material
k_temp(freeWater_domain) = k_freeWater; % assume pure water for cells which consist partly of air and water

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
    air(lowMinOrg_domain) = 1-mineral(lowMinOrg_domain)-organic(lowMinOrg_domain)-water(lowMinOrg_domain);    
    k_temp(lowMinOrg_domain) = (water(lowMinOrg_domain) .* kw.^0.5 + mineral(lowMinOrg_domain) .* km.^0.5 + organic(lowMinOrg_domain).* ko.^0.5 + air(lowMinOrg_domain) .* ka.^0.5).^2 ;
end