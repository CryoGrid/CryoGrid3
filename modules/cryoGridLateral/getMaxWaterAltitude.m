function max_water = getMaxWaterAltitude(PARA)
    max_water = max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater; % L�o : change max->min if you prefer that a relative value of 0, block ponding for everybody.
end
