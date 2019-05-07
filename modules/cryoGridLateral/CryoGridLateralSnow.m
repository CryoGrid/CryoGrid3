function [PARA] = CryoGridLateralSnow( PARA, GRID )
    
    labBarrier();
    % fprintf('\t\t\tsync - calculating snow scaling factors \n');
    
    % calculate scaling factors for all realizations and store in location struct
    PARA.ensemble.snow_scaling = calculateSnowScalingFactors(PARA,GRID);
    PARA.forcing.snow_scaling = PARA.ensemble.snow_scaling(labindex);

    % fprintf('\t\t\tSacling factor of worker %d = %0.3f \n', [labindex, PARA.forcing.snow_scaling] );
end
