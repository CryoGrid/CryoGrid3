function [ wc, GRID ] = removeSediment( GRID, wc )

K_delta = GRID.general.K_delta(GRID.soil.cT_domain);

%assert( GRID.soil.residualOrganic <= 0, 'positive organics to remove' );
%assert( GRID.soil.residualMineral <= 0, 'positive minerals to remove' );

fprintf( '\t\t\tsync - removal of sediment\n' );
fprintf( '\t\t\t\t residual organic:  %3.6e m \n', GRID.soil.residualOrganic );
fprintf( '\t\t\t\t residual mineral:  %3.6e m \n', GRID.soil.residualMineral );

cT_sediment = GRID.soil.cT_mineral + GRID.soil.cT_organic > 1e-9;
firstSedimentCell = find( cT_sediment, 1, 'first' );

sedimentToRemove = K_delta(firstSedimentCell) .* (GRID.soil.cT_mineral(firstSedimentCell) + GRID.soil.cT_organic(firstSedimentCell));

fprintf( '\t\t\t\t sediment to remove:  %3.6e m \n', -sedimentToRemove );



organicToRemove = sedimentToRemove .* double( GRID.soil.residualOrganic<0 ) .* GRID.soil.residualOrganic ./ ...
                   ( double( GRID.soil.residualOrganic<0 ) .* GRID.soil.residualOrganic + double( GRID.soil.residualMineral<0 ) .* GRID.soil.residualMineral ) ;

mineralToRemove = sedimentToRemove .* double( GRID.soil.residualMineral<0 ) .* GRID.soil.residualMineral ./ ... 
                   ( double( GRID.soil.residualOrganic<0 ) .* GRID.soil.residualOrganic + double( GRID.soil.residualMineral<0 ) .* GRID.soil.residualMineral ) ;

removedOrganic = K_delta(firstSedimentCell) .* GRID.soil.cT_organic(firstSedimentCell);
removedMineral = K_delta(firstSedimentCell) .* GRID.soil.cT_mineral(firstSedimentCell);

deltaOrganic = removedOrganic + organicToRemove; % negative if not enough organic in upper cell --> reduce organic in cell below; positive if too much organic removed from upper cell --> increase organic in cell below
deltaMineral = removedMineral + mineralToRemove;

assert( abs( deltaOrganic + deltaMineral ) <= 1e-9, 'residuals organic/mineral do not match' );

GRID.soil.cT_organic(firstSedimentCell)=0;
GRID.soil.cT_mineral(firstSedimentCell)=0;

% adjust fractions in next cell

GRID.soil.cT_organic(firstSedimentCell+1) = max( 0, GRID.soil.cT_organic(firstSedimentCell+1) + deltaOrganic ./ K_delta(firstSedimentCell+1) );
GRID.soil.cT_mineral(firstSedimentCell+1) = max( 0, GRID.soil.cT_mineral(firstSedimentCell+1) + deltaMineral ./ K_delta(firstSedimentCell+1) );

assert( GRID.soil.cT_organic(firstSedimentCell+1)>=0, 'negative organic content' );
assert( GRID.soil.cT_mineral(firstSedimentCell+1)>=0, 'negative mineral content' );

GRID.soil.residualOrganic = min( GRID.soil.residualOrganic + organicToRemove, 0);
GRID.soil.residualMineral = min( GRID.soil.residualMineral + mineralToRemove, 0);

if firstSedimentCell == 1 % no water domain --> remove cells from soil domain
    fprintf( '\t\t\t\t removing cell #%d\n', firstSedimentCell );
    % store water from upper cell
    GRID.soil.water2pool = GRID.soil.water2pool + max( K_delta(firstSedimentCell) .* wc(firstSedimentCell), 0 ); %  take max just to be safe
    
    % adjust air and soil domains and boundaries
    GRID.air.cT_domain(GRID.soil.cT_domain_ub)=1;
    GRID.air.K_domain(GRID.soil.K_domain_ub)=1;
    GRID.air.cT_domain_lb=GRID.air.cT_domain_lb+1;
    GRID.air.K_domain_lb=GRID.air.K_domain_lb+1;
    GRID.soil.cT_domain(GRID.soil.cT_domain_ub)=0;
    GRID.soil.K_domain(GRID.soil.K_domain_ub)=0;
    GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub+1;
    GRID.soil.K_domain_ub=GRID.soil.K_domain_ub+1;
    GRID.soil.soilGrid(firstSedimentCell)=[];
    
    wc(firstSedimentCell)=[];
    GRID.soil.cT_water(firstSedimentCell)=[];
    
    GRID.soil.cT_organic(firstSedimentCell)=[];
    GRID.soil.cT_natPor(firstSedimentCell)=[];
    GRID.soil.cT_actPor(firstSedimentCell)=[];
    GRID.soil.cT_mineral(firstSedimentCell)=[];
    GRID.soil.cT_soilType(firstSedimentCell)=[];
    
    GRID.soil.excessGroundIce(firstSedimentCell)=[];
    
else % water domain present
    fprintf( '\t\t\t\t water above soil\n' )
    fprintf( '\t\t\t\t removing cell #%d\n', firstSedimentCell );

    % adjust actual porosity and soil type ( natural porosity unchanged )
    GRID.soil.cT_actPor(firstSedimentCell) = 1;
    GRID.soil.cT_soilType(firstSedimentCell) = 3; % set to water body to have correct field capacity
    
    % move water down / air up
    totalWater = sum( wc(1:firstSedimentCell) .* K_delta(1:firstSedimentCell) );
    k = firstSedimentCell;
    while k>0
        %assert( totalWater >= 0, sprintf( 'negative water to pond in cell #%d', k ) );
        
        tempWater= max( 0, min( totalWater, K_delta(k) .* GRID.soil.cT_actPor(k) ) );
        wc(k)=tempWater./K_delta(k);
        totalWater=totalWater-tempWater;
        k=k-1;
    end
    assert( abs(totalWater) < 1e-6, sprintf( '\t\t\t\t\t\t water to pond after movin water down: totalWater=%3.6e\n', totalWater ) );

    % check if uppermost cell of soil domain contains only air --> remove
    while (GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)+wc(1)<1e-6)
        fprintf( '\t\t\t\t removing air cell\n' )

        % adjust air and soil domains and boundaries
        GRID.air.cT_domain(GRID.soil.cT_domain_ub)=1;
        GRID.air.K_domain(GRID.soil.K_domain_ub)=1;
        GRID.air.cT_domain_lb=GRID.air.cT_domain_lb+1;
        GRID.air.K_domain_lb=GRID.air.K_domain_lb+1;
        GRID.soil.cT_domain(GRID.soil.cT_domain_ub)=0;
        GRID.soil.K_domain(GRID.soil.K_domain_ub)=0;
        GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub+1;
        GRID.soil.K_domain_ub=GRID.soil.K_domain_ub+1;
        GRID.soil.soilGrid(1)=[];

        wc(1)=[];
        GRID.soil.cT_water(1)=[];
        
        GRID.soil.cT_organic(1)=[];
        GRID.soil.cT_natPor(1)=[];
        GRID.soil.cT_actPor(1)=[];
        GRID.soil.cT_mineral(1)=[];
        GRID.soil.cT_soilType(1)=[];

        GRID.soil.excessGroundIce(1)=[];

    end
        
    % just to be sure
    GRID.soil.cT_water = wc;
    
end

end
