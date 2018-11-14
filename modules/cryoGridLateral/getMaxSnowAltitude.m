function max_snow = getMaxSnowAltitude(PARA)
    % max_snow = max(PARA.ensemble.altitude) + PARA.snow.relative_maxSnow;
    max_snow = PARA.ensemble.altitude(labindex) + PARA.snow.relative_maxSnow; % Leo's veresion to check pb when too much snow
end
