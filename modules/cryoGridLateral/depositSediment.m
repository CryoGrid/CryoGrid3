function [ wc, GRID ] = depositSediment( GRID, wc )

K_delta = GRID.general.K_delta(GRID.soil.cT_domain);

%assert( GRID.soil.residualOrganic >= 0, 'negative organics to deposit' );
%assert( GRID.soil.residualMineral >= 0, 'negative minerals to deposit' );

fprintf( '\t\t\tsync - deposition of sediment\n' );
fprintf( '\t\t\t\t residual organic:  %3.6e m \n', GRID.soil.residualOrganic );
fprintf( '\t\t\t\t residual mineral:  %3.6e m \n', GRID.soil.residualMineral );

cT_sediment = GRID.soil.cT_mineral + GRID.soil.cT_organic > 1e-9;
firstSedimentCell = find( cT_sediment, 1, 'first' );

% determine amount to deposit
sedimentToDeposit = K_delta(firstSedimentCell) .* ( 1 - GRID.soil.cT_natPor(firstSedimentCell) - GRID.soil.cT_mineral(firstSedimentCell) - GRID.soil.cT_organic(firstSedimentCell) );
depositionType = 1;
if abs(sedimentToDeposit)<1e-9
    sedimentToDeposit = K_delta(firstSedimentCell) .* ( 1 - GRID.soil.cT_natPor(firstSedimentCell) );
    depositionType = 2;
end

depositedFractionOrganic = double( GRID.soil.residualOrganic>0 ) .* GRID.soil.residualOrganic ./ ...
                            ( double( GRID.soil.residualOrganic>0 ) .* GRID.soil.residualOrganic + double( GRID.soil.residualMineral>0 ) .* GRID.soil.residualMineral );

depositedFractionMineral = double( GRID.soil.residualMineral>0 ) .* GRID.soil.residualMineral ./ ...
                            ( double( GRID.soil.residualOrganic>0 ) .* GRID.soil.residualOrganic + double( GRID.soil.residualMineral>0 ) .* GRID.soil.residualMineral );

                        

% case distinction: filling uppermost sediment cell versus creating new cell
if depositionType == 1  % filling upppermost sediment cell
    fprintf( 'filling first sediment cell\n' );
    

    
    GRID.soil.cT_organic(firstSedimentCell) = GRID.soil.cT_organic(firstSedimentCell) + sedimentToDeposit .* depositedFractionOrganic ./ K_delta(firstSedimentCell);
    GRID.soil.cT_mineral(firstSedimentCell) = GRID.soil.cT_mineral(firstSedimentCell) + sedimentToDeposit .* depositedFractionMineral ./ K_delta(firstSedimentCell);
    GRID.soil.cT_actPor(firstSedimentCell) = 1 - GRID.soil.cT_organic(firstSedimentCell) - GRID.soil.cT_mineral(firstSedimentCell);
    
    GRID.soil.residualOrganic = max( GRID.soil.residualOrganic - sedimentToDeposit .* depositedFractionOrganic, 0);
    GRID.soil.residualMineral = max( GRID.soil.residualMineral - sedimentToDeposit .* depositedFractionMineral, 0);
    
    excessWater = max( 0, K_delta(firstSedimentCell) .* ( wc(firstSedimentCell) - GRID.soil.cT_actPor(firstSedimentCell) ) );
    GRID.soil.water2pool = GRID.soil.water2pool + excessWater;
    wc(firstSedimentCell) = wc(firstSedimentCell) - excessWater ./ K_delta(firstSedimentCell);
    GRID.soil.cT_water(firstSedimentCell) = wc(firstSedimentCell);
    
    assert( GRID.soil.cT_mineral(firstSedimentCell)+GRID.soil.cT_organic(firstSedimentCell)+wc(firstSedimentCell) <= 1, 'first sediment cell exceeded' );
    
    
elseif depositionType == 2
    fprintf( 'creating new sediment cell - ' );
    % case distinction depending on water domain present or not
    if firstSedimentCell==1    % no water body --> create new soil cell
        fprintf( 'within air domain\n' );
        
        % create new soil cell / change GRID domains
        GRID.soil.cT_domain(GRID.air.cT_domain_lb)=1;
        GRID.soil.K_domain(GRID.air.K_domain_lb)=1;
        GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub-1;
        GRID.soil.K_domain_ub=GRID.soil.K_domain_ub-1;
        GRID.air.cT_domain(GRID.air.cT_domain_lb)=0;
        GRID.air.K_domain(GRID.air.K_domain_lb)=0;
        GRID.air.cT_domain_lb=GRID.air.cT_domain_lb-1;
        GRID.air.K_domain_lb=GRID.air.K_domain_lb-1;
        
        % fill new soil cell
        cellSize = K_delta(firstSedimentCell);
        GRID.soil.cT_natPor =  [ GRID.soil.cT_natPor(firstSedimentCell); GRID.soil.cT_natPor ];    % take natPor of cell below
        
        %potentialSedimentNewCell = K_delta(1) .* ( 1 - GRID.soil.cT_natPor(1) );
        
        %addedSediment = min( GRID.soil.residualOrganic+GRID.soil.residualMineral,  potentialSedimentNewCell );
        
        GRID.soil.cT_organic =  [ depositedFractionOrganic .* sedimentToDeposit ./ cellSize ; GRID.soil.cT_organic ];
        GRID.soil.cT_mineral =  [ depositedFractionMineral .* sedimentToDeposit ./ cellSize ; GRID.soil.cT_mineral ];
        GRID.soil.cT_actPor =   [ 1-GRID.soil.cT_organic(1)-GRID.soil.cT_mineral(1) ; GRID.soil.cT_actPor ];
        
        GRID.soil.residualOrganic = max( GRID.soil.residualOrganic - sedimentToDeposit .* depositedFractionOrganic, 0);
        GRID.soil.residualMineral = max( GRID.soil.residualMineral - sedimentToDeposit .* depositedFractionMineral, 0);
        
        % update remaining soil fields
        wc = [ 0 ; wc ];
        GRID.soil.cT_water = [ 0; GRID.soil.cT_water ];
        GRID.soil.cT_soilType = [ GRID.soil.cT_soilType(firstSedimentCell); GRID.soil.cT_soilType];
        GRID.soil.excessGroundIce = [ 0 ; GRID.soil.excessGroundIce ];
        
        % update GRID spacings
        GRID.general.K_grid(GRID.soil.cT_domain_ub) = GRID.general.K_grid(GRID.soil.cT_domain_ub+1)-cellSize;
        GRID.general.K_grid(GRID.air.cT_domain) = [GRID.general.K_grid(GRID.air.cT_domain_lb)+(-GRID.snow.snowCellSize)*(GRID.air.cT_domain_lb-1):GRID.snow.snowCellSize:GRID.general.K_grid(GRID.air.cT_domain_lb)]';
        GRID.general.cT_grid = ( GRID.general.K_grid(1:end-1)+ GRID.general.K_grid(2:end))/2;
        GRID.general.cT_delta = (- GRID.general.cT_grid(1:end-1,1)+ GRID.general.cT_grid(2:end,1));
        GRID.general.K_delta = (- GRID.general.K_grid(1:end-1,1)+ GRID.general.K_grid(2:end,1));
        GRID.soil.soilGrid = [ GRID.general.K_grid(GRID.soil.cT_domain_ub) ; GRID.soil.soilGrid ];
        
    elseif firstSedimentCell>1      % water body --> deposit in water cell
        fprintf( 'within water domain\n' );
        
        lastWaterCell = find( GRID.soil.cT_organic + GRID.soil.cT_mineral <= 1e-9, 1, 'last' );
        
        assert( lastWaterCell == firstSedimentCell-1, 'water / sediment boundray not clear' );
        
        fprintf( 'depositing sediment in water cell # %d\n', lastWaterCell );

        GRID.soil.cT_organic(lastWaterCell) =  depositedFractionOrganic .* sedimentToDeposit ./ K_delta(lastWaterCell);
        GRID.soil.cT_mineral(lastWaterCell) =  depositedFractionMineral .* sedimentToDeposit ./ K_delta(lastWaterCell);
        GRID.soil.cT_actPor(lastWaterCell) = 1 - GRID.soil.cT_organic(lastWaterCell) - GRID.soil.cT_mineral(lastWaterCell);
        
        % inherit soil type from uppermost sediment cell
        GRID.soil.cT_soilType(lastWaterCell) = GRID.soil.cT_soilType(firstSedimentCell);
        
        % update residual storages
        GRID.soil.residualOrganic = max( GRID.soil.residualOrganic - sedimentToDeposit .* depositedFractionOrganic, 0);
        GRID.soil.residualMineral = max( GRID.soil.residualMineral - sedimentToDeposit .* depositedFractionMineral, 0);

        excessWater = max( 0, K_delta(lastWaterCell) .* ( wc(lastWaterCell) - GRID.soil.cT_actPor(lastWaterCell) ) );
        GRID.soil.water2pool = GRID.soil.water2pool + excessWater;
        
        wc(lastWaterCell) = wc(lastWaterCell) - excessWater ./ K_delta(lastWaterCell);
        GRID.soil.cT_water(lastWaterCell) = wc(lastWaterCell);
        
        assert( GRID.soil.cT_mineral(lastWaterCell)+GRID.soil.cT_organic(lastWaterCell)+wc(lastWaterCell) <= 1, 'last water cell exceeded' );
        
    end
end