function max_snow = getMaxSnowAltitude(PARA)
    max_snow = PARA.ensemble.altitude(labindex) + PARA.snow.relative_maxSnow;
end
