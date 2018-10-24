function max_snow = getMaxSnowAltitude(PARA)
<<<<<<< HEAD
    max_snow = PARA.ensemble.altitude(labindex) + PARA.snow.relative_maxSnow;
=======
    max_snow = max(PARA.ensemble.altitude) + PARA.snow.relative_maxSnow;
>>>>>>> origin/xice_mpi_polygon_TC
end
