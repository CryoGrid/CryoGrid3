function max_water = getMaxWaterAltitude(PARA)
    max_water = max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater;
end
