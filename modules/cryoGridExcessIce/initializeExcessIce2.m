function [GRID,PARA] = initializeExcessIce2(GRID,PARA) 

GRID.soil.excessGroundIce = GRID.soil.cT_water>GRID.soil.cT_natPor;