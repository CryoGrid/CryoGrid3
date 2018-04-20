function max_water = getMaxWaterAltitude(PARA)
    max_water = min( PARA.ensemble.soil_altitude ) + PARA.soil.relative_maxWater;
end
