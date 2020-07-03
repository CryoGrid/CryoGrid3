function snow_change = calculateLateralSnowFluxes( mobile_snow, PARA )

	index = labindex;
    % calculate mobile snow volume
    mobile_snow_volume = sum ( mobile_snow .* PARA.ensemble.weight .* (PARA.ensemble.terrain_index_snow>0) );
    % calculate individual change (deposit or removal) dependent on terrain index
    if PARA.ensemble.terrain_index_snow(index)>0 %removal
        snow_change = -mobile_snow(index);
    elseif PARA.ensemble.terrain_index_snow(index)<0 % deposit
        snow_change = -mobile_snow_volume .* PARA.ensemble.terrain_index_snow(index);
    else % terrain index = 0 or NaN
        snow_change = 0;
    end
    
end
