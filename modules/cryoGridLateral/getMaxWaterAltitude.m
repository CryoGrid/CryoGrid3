function max_water = getMaxWaterAltitude(PARA)
%    max_water = max( PARA.ensemble.soil_altitude ) + PARA.soil.relative_maxWater;
%tsvd IS use minimum instead of maximum for IS 
    %max_water = min( PARA.ensemble.soil_altitude ) + PARA.soil.relative_maxWater; % assume that ponding only happens on lower tundra tile - define relative_maxWater w.r.t. this lower level (i.e. use min instead of max in function)
    max_water = PARA.ensemble.soil_altitude(end) + PARA.soil.relative_maxWater; % set to ponding value of last tundra tile
end
