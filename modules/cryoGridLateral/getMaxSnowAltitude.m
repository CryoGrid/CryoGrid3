function max_snow = getMaxSnowAltitude(PARA)
    max_snow = max( PARA.ensemble.altitude ) + PARA.snow.relative_maxSnow;
end
