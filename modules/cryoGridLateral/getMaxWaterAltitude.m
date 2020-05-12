function max_water = getMaxWaterAltitude(PARA)
    max_water = max( PARA.ensemble.soil_altitude ) + PARA.soil.relative_maxWater;
end
