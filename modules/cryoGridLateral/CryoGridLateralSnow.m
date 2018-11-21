function [T, GRID, BALANCE, TEMPORARY] = CryoGridLateralSnow( PARA, GRID, BALANCE, TEMPORARY, FORCING, T)
    
    labBarrier();
    fprintf('\t\t\tsync - calculating snow scaling factors \n');
    
    factors = calculateSnowScalingFactors(PARA,GRID);
    PARA.forcing.snow_scaling = factors(labindex);
    %
    %BALANCE.water.dr_lateralSnow = BALANCE.water.dr_lateralSnow + my_snow_change*1000 ;
    %TEMPORARY.snow_flux_lateral = TEMPORARY.snow_flux_lateral + my_snow_change;
    
    fprintf('\t\t\tSnow flux to worker %d = %f mm SWE \n', [labindex, my_snow_change*1000] );
end
