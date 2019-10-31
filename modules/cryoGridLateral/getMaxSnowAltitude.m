function max_snow = getMaxSnowAltitude(PARA)
    isInMyComp = ( conncomp(graph(PARA.ensemble.adjacency_snow)) == labindex );  
    max_snow = max(PARA.ensemble.altitude( isInMyComp )) + PARA.snow.relative_maxSnow;
end
