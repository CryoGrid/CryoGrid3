function c_temp = capacityUnfrozen(wc,GRID,PARA)

c_w = PARA.constants.c_w; % 4.2*10^6; %[J/m�K]
c_o = PARA.constants.c_o; % 2.5*10^6; %[J/m�K]
c_m = PARA.constants.c_m; % 2*10^6; %[J/m�K]

c_temp = (GRID.soil.cT_mineral+GRID.soil.cT_organic>1e-6) .* (GRID.soil.cT_mineral.*c_m + GRID.soil.cT_organic.*c_o + wc.*c_w) + ...
         (GRID.soil.cT_mineral+GRID.soil.cT_organic<=1e-6) .* c_w ; % assume pure water for cells which consist partly of air and water
