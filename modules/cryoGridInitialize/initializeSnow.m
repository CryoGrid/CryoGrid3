function GRID = initializeSnow(GRID)

GRID.snow.Snow_i=zeros(size(GRID.air.cT_domain));
GRID.snow.Snow_w=zeros(size(GRID.air.cT_domain));
GRID.snow.Snow_a=zeros(size(GRID.air.cT_domain));
GRID.snow.SWEinitial=0;