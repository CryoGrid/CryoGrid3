function [ flag,finalThick ] = excessGroundIceCheckThickness(GRID, threshold)
% When Xice + lateral heat are on, this function checks if cells are thick
% enough to support shrinking when the Xice module cannot route the water
% up because of overlying frozen cells.

cellThick = GRID.general.K_delta(GRID.soil.cT_domain);

finalThick = cellThick .* (GRID.soil.cT_mineral + GRID.soil.cT_organic) ./ (1-GRID.soil.cT_natPor); % Caltuclate thickness of cell if skrinked by excess ice melt

if sum(finalThick < threshold) > 0; % Check if some cells are thinner than threshold
    flag = 0;
else
    flag = 1;
    
end

end