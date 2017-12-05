function k_temp = capacityUnfrozen(wc,GRID, PARA)

ka = PARA.constants.k_a; %0.025;       %air [Hillel(1982)]
kw = PARA.constants.k_w; %0.57;        %water [Hillel(1982)]
ko = PARA.constants.k_o; %0.25;        %organic [Hillel(1982)]
km = PARA.constants.k_m; %soil.kh_bedrock;     %mineral 

air=1-wc-GRID.soil.cT_mineral-GRID.soil.cT_organic;

k_temp = (GRID.soil.cT_mineral+GRID.soil.cT_organic>1e-6) .* (wc.* kw.^0.5 + GRID.soil.cT_mineral.* km.^0.5 + GRID.soil.cT_organic.* ko.^0.5 + air.* ka.^0.5).^2 + ...
         (GRID.soil.cT_mineral+GRID.soil.cT_organic<=1e-6) .* kw; % assume pure water for cells which consist partly of air and water


