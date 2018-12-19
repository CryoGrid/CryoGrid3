function [ wc, GRID ] = CryoGridLateralErosion( PARA, GRID, wc, T )

labBarrier();
% check preconditions (are there two connected realizations which
% between which erosion occurs
precondition_sedimentExchange = checkPreconditionWaterExchange( T, GRID ); % identical to water exchange: unfrozen upp


if precondition_sedimentExchange
    
    fprintf('\t\t\tsync - exchanging sediment (lateral erosion)\n');
    
    sediment_change_tot = zeros(1,numlabs); % in [m] w.r.t. current realization (index)
    sediment_change_diff = zeros(1,numlabs);
    sediment_change_adv = zeros(1,numlabs);
    sediment_change_o = zeros(1,numlabs);
    sediment_change_m = zeros(1,numlabs);
    
    
    PACKAGE_sedimentExchange.K_delta = GRID.general.K_delta( GRID.soil.cT_domain );
    PACKAGE_sedimentExchange.cT_mineral = GRID.soil.cT_mineral;
    PACKAGE_sedimentExchange.cT_organic = GRID.soil.cT_organic;
    PACKAGE_sedimentExchange.erosion_condition = T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub);    %identical to water exchange
    
    for j=1:numlabs
        if j~=labindex
            labSend( PACKAGE_sedimentExchange, j, 2);
        end
    end
    for j=1:numlabs
        if j~=labindex
            PACKAGE_sedimentExchange_j = labReceive(j, 2);
            [   sediment_change_tot(j),...
                sediment_change_diff(j),...
                sediment_change_adv(j),...
                sediment_change_o(j),...
                sediment_change_m(j),...
                GRID                      ] = calculateLateralSedimentFluxes( PARA, GRID, PACKAGE_sedimentExchange_j, T, j );
            
        end
    end
    
    fprintf('\t\t\tsync - nominal total sediment flux to worker %d = %3.6e m \n', [labindex, sum(sediment_change_tot) ] );
    fprintf('\t\t\tsync - diffusive sediment flux to worker %d = %3.6e m \n', [labindex, sum(sediment_change_diff) ] );
    fprintf('\t\t\tsync - advective sediment flux to worker %d = %3.6e m \n', [labindex, sum(sediment_change_adv) ] );
    fprintf('\t\t\tsync - organic sediment flux to worker %d = %3.6e m \n', [labindex, sum(sediment_change_o) ] );
    fprintf('\t\t\tsync - mineral sediment flux to worker %d = %3.6e m \n', [labindex, sum(sediment_change_m) ] );
    
    GRID.soil.residualOrganic = GRID.soil.residualOrganic + sum( sediment_change_o );
    GRID.soil.residualMineral = GRID.soil.residualMineral + sum( sediment_change_m );
    GRID.soil.residualSediment = GRID.soil.residualSediment + sum( sediment_change_tot );
    
    if T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub)
        % determine uppermost cell containing sediment
        cT_sediment = GRID.soil.cT_mineral + GRID.soil.cT_organic > 1e-9;
        firstSedimentCell = find( cT_sediment, 1, 'first' );
        
        K_delta = GRID.general.K_delta(GRID.soil.cT_domain);
        
        % determine thresholds for removal or creation of cell
        thresholdRemoveCell = K_delta(firstSedimentCell) .* (GRID.soil.cT_mineral(firstSedimentCell) + GRID.soil.cT_organic(firstSedimentCell));
        thresholdNewCell = K_delta(firstSedimentCell) .* (1 - GRID.soil.cT_natPor(firstSedimentCell));
        applyFluxes =  GRID.soil.residualSediment >= thresholdNewCell || GRID.soil.residualSediment <= -thresholdRemoveCell ;
        
        % apply fluxes if needed (for now only this is done by each realization individually
        if applyFluxes
            fprintf( '\t\t\tsync - applying sediment fluxes to worker %d\n ...', labindex );          
            
            if GRID.soil.residualSediment >= thresholdNewCell
                % ------  deposition of sediment------------------------
                [ wc, GRID ] = depositSediment( PARA, GRID, wc );
                
                assert( legnth(wc)==sum(GRID.soil.cT_domain), 'wc length wrong after sediment deposition');
                
            elseif GRID.soil.residualSediment <= -thresholdRemoveCell
                % ------ removal of sediment -----------------
                [ wc, GRID ] = removeSediment( GRID, wc );
                
                assert( legnth(wc)==sum(GRID.soil.cT_domain), 'wc length wrong after sediment removal');

            else
                error('lateral erosion - sediment flux wrongly applied');
            end
            
            % track total amount (volume) of minerals and organics in order to
            % confirm conservation of mass --> BALANCE.mass ?
            
            GRID.soil.residualSediment = 0;
            assert( GRID.soil.residualOrganic == 0, 'organics not zero after application of fluxes');
            assert( GRID.soil.residualMineral == 0, 'minerals not zero after application of fluxes');
            
            % update water body domain and LUT
            GRID = updateGRID_erosion( PARA, GRID );
        end       
    end
end

end