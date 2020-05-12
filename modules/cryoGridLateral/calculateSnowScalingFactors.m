function factors = calculateSnowScalingFactors(PARA,GRID)
    weights = PARA.ensemble.weight; % relative areal weights of all tiles, arbitratry unit
    immobileSnowAltitudes = PARA.ensemble.altitude + PARA.ensemble.immobile_snow_height;
    surfaceAltitudes = PARA.ensemble.surface_altitude;
    
    assert( isrow(immobileSnowAltitudes) && isrow(surfaceAltitudes), 'vector alignment is not correct');
    
    delta = GRID.snow.snowCellSize;
    
    % determine workers in same connected component
    bins = conncomp(graph(PARA.ensemble.adjacency_snow));
    isInMyComp = ( double(bins) == bins(labindex) );  
    % rescale snow input for each component of the current worker    
    factors = zeros(1, numlabs);      
    for i=find( isInMyComp ) % iterate over the component of the current worker

        % determine number of lower-lying tiles
        lowerLyingInComp = ( surfaceAltitudes + delta < surfaceAltitudes(i) ) & ( isInMyComp );

        if sum(double(lowerLyingInComp))>0 && surfaceAltitudes(i)>=immobileSnowAltitudes(i)   % in case there is mobile snow and lower-lying tiles --> distribute
            factors(lowerLyingInComp) = factors(lowerLyingInComp) + weights(i) ./ sum( weights(lowerLyingInComp) );
        else % otherwise don't distribute
            factors(i) = factors(i) + 1.0;
        end
        
    end
    
    if sum( weights(isInMyComp) ) ~= sum( factors(isInMyComp) .* weights(isInMyComp) )
        warning('Snow exchange is not mass-conserving!');
    end

end