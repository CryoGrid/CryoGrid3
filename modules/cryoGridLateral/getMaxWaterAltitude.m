function max_water = getMaxWaterAltitude(PARA)
    max_water = min( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater; % Léo : max->min so that a relative value of 0, block ponding for everybody.
end
