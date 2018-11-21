function factors = calculateSnowScalingFactors()
    weights = PARA.ensemlbe.weight; % relative areal weights of all tiles, arbitratry unit
    maxSnowAltitude = PARA.location.absolute_maxSnow_altitude;
    immobileSnowAltitudes = PARA.ensemble.altitude + PARA.ensemble.immobile_snow_height;
    surfaceAltitudes = PARA.ensemble.surface_altitude;
    
    delta = GRID.snow.snowCellSize;
    
    factors = zeros(1, numlabs);      
    for i=1:numlabs
        
        if surfaceAltitudes(i)<maxSnowAltitude
            
            % determine number of lower-lying tiles
            lowerLying = surfaceAltitudes + delta < surfaceAltitudes(i);
            
            if sum(double(lowerLying))>0 && surfaceAltitudes(i)>=immobileSnowAltitudes(i)   % in case there is mobile snow and lower-lying tiles --> distribute
                factors(lowerLying) = factors(lowerLying) + weights(i) ./ sum( weights(lowerLying) );
            else % otherwise don't distribute
                factors(i) = factors(i) + 1.0;
            end
        
        else
            % let the factor zero, store the removed excess snow
        end
    end
    
    if sum( weights ) ~= sum( factors .* weights )
        warning('Snow exchange is not mass-conserving!');
    end


end