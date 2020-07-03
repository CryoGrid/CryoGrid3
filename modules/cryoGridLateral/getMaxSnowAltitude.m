function max_snow = getMaxSnowAltitude(PARA)
% %tsvd IS  % allow for twice snow accumulation at embankment shoulder    
%     if(PARA.IS.tile_type=='shoulder_ag')
%         max_snow = max(PARA.ensemble.altitude) + 2*PARA.snow.relative_maxSnow; 
%     else
        max_snow = max(PARA.ensemble.altitude) + PARA.snow.relative_maxSnow;
%    end
    
end
