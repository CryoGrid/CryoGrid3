function max_snow = getMaxSnowAltitude(PARA)
    bins = conncomp(graph(PARA.ensemble.adjacency_snow));
    isInMyComp = ( double(bins) == bins(labindex) );
    max_snow = max(PARA.ensemble.altitude( isInMyComp )) + PARA.snow.relative_maxSnow;
end
