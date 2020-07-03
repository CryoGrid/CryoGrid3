function [wc, GRID, BALANCE] = CryoGridLateralWater( PARA, GRID, BALANCE, T, wc)

    % labBarrier();
    %tsvd IS  allocate variables here
    water_fluxes        = zeros(numlabs,numlabs); % in m of height change
    water_fluxes_worker = zeros(numlabs,numlabs); % account for available water
    water_fluxes_gather = zeros(numlabs,numlabs,numlabs); 
    boundary_water_flux = 0.; 

%%%    disp(['check precondH2O!!! for worker...',num2str(labindex)])
    % check preconditions
    precondition_waterExchange = checkPreconditionWaterExchange( T, GRID );

        % WRAPPER
    fprintf('\t\t\tsync - exchanging water\n');
    labBarrier();
    % calculate lateral water fluxes
    %tsvd IS        water_fluxes= zeros(1,numlabs); % in m of height change  - use vector instead of matrix
    PACKAGE_waterExchange.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
    PACKAGE_waterExchange.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
    PACKAGE_waterExchange.infiltration_condition = T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub);

    j_prevWorker = mod(labindex - 2, numlabs) + 1; j_nextWorker = mod(labindex, numlabs) + 1; 
    j_adjacent = [j_prevWorker,j_nextWorker];
    PACKAGE_waterExchange_prev = labSendReceive(j_nextWorker,j_prevWorker,PACKAGE_waterExchange,2); PACKAGE_waterExchange_next = labSendReceive(j_prevWorker,j_nextWorker,PACKAGE_waterExchange,2);
    PACKAGE_waterExchange_adjacent = [PACKAGE_waterExchange_prev,PACKAGE_waterExchange_next];

    if (precondition_waterExchange) 
        for jj=1:2 % lateral water exchange with previous and next tile
            j = j_adjacent(jj); % first exchange with previous tile, then exchange with next tile     
            PACKAGE_waterExchange_j = PACKAGE_waterExchange_adjacent(jj);
    %         [F_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/s ]  
            water_fluxes = calculateLateralWaterDarcyFluxes( T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes);  % matrix containing all fluxes in [m/s] scaled to row index
            if(numlabs==2) % if only 2 tiles, only 1 loop is needed
                break
            end
        end
    end
    %%% check if water_fluxes is updated correctly by second loop!  AND  prevent fluxes for outmost tiles...
    
%     j_prevWorker = mod(labindex - 2, numlabs) + 1;
%     j_nextWorker = mod(labindex, numlabs) + 1;
%     PACKAGE_waterExchange_j = labSendReceive(j_nextWorker,j_prevWorker,PACKAGE_waterExchange,2)
%     if (precondition_waterExchange && labindex>1) % worker 1 not for lateral exchange (no lower worker connected)
%         j=labindex-1; % number of worker for lateral exchange with current worker(i.e. exchange with the next lower worker) 
%         water_fluxes = calculateLateralWaterDarcyFluxes( T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes);  % matrix containing all fluxes in [m/s] scaled to row index
%     end  
        
    fprintf('\t\t\tsync - exchanging water   end calc lat water darcy fluxes\n');  %tttt
        
    if(T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub)) %tsvd IS   allow boundary fluxes for individual tiles (no assumption that workers are connected)
        % Calculate possible boundary fluxes
        [ boundary_water_flux ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T);
        fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
    end
    % Check for water availability and set real water fluxes
    [ water_fluxes_worker, boundary_water_flux ] = calculateLateralWaterAvailable( PARA,GRID, wc, water_fluxes,boundary_water_flux );
    fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)

        water_fluxes_gather(:,:,labindex)=water_fluxes_worker; % allows to account for update water fluxes when not water available! This information needs calculation on each worker separately...
        % Send real fluxes all around
        for j=1:numlabs  %ccc restrict to pairwise!?
            if j~=labindex
                labSend( water_fluxes_worker, j, 5);
            end
        end

        for j=1:numlabs
            if j~=labindex
                water_fluxes_worker_j = labReceive(j, 5);
                water_fluxes_gather(:,:,j)=water_fluxes_worker_j;
            end
        end

        water_fluxes=nansum(water_fluxes_gather,3);
        waterflux=nansum(water_fluxes(labindex,:))+boundary_water_flux;

        % waterflux=nansum(water_fluxes_worker(labindex,:))+boundary_water_flux;  % fix Jan...  
                 
%    water_fluxes_gather(:,:,labindex) = water_fluxes_worker;     
    % Send real fluxes all around       
% 	if(labindex<numlabs)
% 		labSend( water_fluxes_worker, labindex+1, 5);
% 		fprintf('\t\t\t ccc labSend ok! \n');      
% 	end
% 
%     if(labindex>1)
%         water_fluxes_worker_j = labReceive(labindex-1, 5);
%         water_fluxes_gather(:,:,j)=water_fluxes_worker_j;
%     end

%     water_fluxes_j = labSendReceive(j_nextWorker,j_prevWorker,water_fluxes_worker,5)
%     if labindex>1
%         water_fluxes_gather(:,:,j_prevWorker)=water_fluxes_j; %tsvd IS  only flux for next lower worker is non-zero!  ...ccc
%     end
    disp('water fluxes: ')
%%% ccc    water_fluxes=nansum(water_fluxes_gather,3) 
   
% outcomment....
%    water_fluxes=water_fluxes_worker;
%    waterflux=nansum(water_fluxes(labindex,:))+boundary_water_flux;
        
    % apply lateral water flux directly (as bulk subsurface flux)
    [wc, excess_water, lacking_water] = bucketScheme(T, wc, zeros( size(wc) ), GRID, PARA, waterflux);
    try
        assert( lacking_water < 1e-9, 'CryoGrid3 - lateral exchange - lacking water>0');    % there should be no lacking water as this was checked for
    catch
        fprintf('\t\t\tLacking water = %3.2e m\n', lacking_water );
    end

    % Store and display
    BALANCE.water.dr_water_fluxes_out=BALANCE.water.dr_water_fluxes_out+water_fluxes./(PARA.technical.syncTimeStep*24*3600); % here we decide in which units we want it. Now m/s
    BALANCE.water.dr_lateralWater = BALANCE.water.dr_lateralWater + (waterflux-excess_water)*1000; % Excess water is removed so that we only keep the net water modification implied by the lateral fluxes
    fprintf('\t\t\t Net wc change :\t%3.2e m\n',waterflux-excess_water)
    if excess_water>1e-9
        GRID.soil.water2pool= GRID.soil.water2pool + excess_water;
        fprintf('\t\t\tExcess water :\t%3.2e m\n',excess_water)
        BALANCE.water.dr_lateralExcess=BALANCE.water.dr_lateralExcess + excess_water*1000;            % Added by Leo to have the lateral fluxes in BALANCE
    end
    if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1
        BALANCE.water.dr_DarcyReservoir = BALANCE.water.dr_DarcyReservoir + boundary_water_flux*1000;
        fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
    end
   % end % end if condition
   % end % end precondition

end

% %%% Original:  %%%%%%%%%%%%%%%%%
% 
% function [wc, GRID, BALANCE] = CryoGridLateralWater( PARA, GRID, BALANCE, T, wc)
% 
%     labBarrier();
%     % check preconditions
%     precondition_waterExchange = checkPreconditionWaterExchange( T, GRID );
%     if precondition_waterExchange
%         % WRAPPER
%         fprintf('\t\t\tsync - exchanging water\n');
%         % calculate lateral water fluxes
%         water_fluxes= zeros(numlabs,numlabs); % in m of height change
%         PACKAGE_waterExchange.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
%         PACKAGE_waterExchange.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
%         PACKAGE_waterExchange.infiltration_condition = T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub);
%         PACKAGE_waterExchange.soil_altitude = PARA.ensemble.soil_altitude(labindex);
% 
%         for j=1:numlabs
%             if j~=labindex
%                 labSend( PACKAGE_waterExchange, j, 2);
%             end
%         end
%         for j=1:numlabs
%             if j~=labindex
%                 PACKAGE_waterExchange_j = labReceive(j, 2);
%                 water_fluxes = calculateLateralWaterDarcyFluxes( T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes);  % matrix containing all fluxes in [m/s] scaled to row index
%             end
%         end
% 
%         % Calculate possible boundary fluxes
%         [ boundary_water_flux ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T);
%         %fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
% 
% 
%         % Check for water availability and set real water fluxes
%         [ water_fluxes_worker, boundary_water_flux ] = calculateLateralWaterAvailable( PARA,GRID, wc, water_fluxes,boundary_water_flux );
%         water_fluxes_gather=zeros(numlabs,numlabs,numlabs);
%         water_fluxes_gather(:,:,labindex)=water_fluxes_worker;
%         %fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
% 
% 
%         % Send real fluxes all around
%         for j=1:numlabs
%             if j~=labindex
%                 labSend( water_fluxes_worker, j, 5);
%             end
%         end
% 
%         for j=1:numlabs
%             if j~=labindex
%                 water_fluxes_worker_j = labReceive(j, 5);
%                 water_fluxes_gather(:,:,j)=water_fluxes_worker_j;
%             end
%         end
% 
%         water_fluxes=nansum(water_fluxes_gather,3);
%         waterflux=nansum(water_fluxes(labindex,:))+boundary_water_flux;
% 
%         % apply lateral water flux directly (as bulk subsurface flux)
%         [wc, excess_water, lacking_water] = bucketScheme(T, wc, zeros( size(wc) ), GRID, PARA, waterflux);
%         try
%             assert( lacking_water < 1e-9, 'CryoGrid3 - lateral exchange - lacking water>0');    % there should be no lacking water as this was checked for
%         catch
%             fprintf('\t\t\tLacking water = %3.2e m\n', lacking_water );
%         end
% 
%         % Store and display
%         BALANCE.water.dr_water_fluxes_out=BALANCE.water.dr_water_fluxes_out+water_fluxes./(PARA.technical.syncTimeStep*24*3600); % here we decide in which units we want it. Now m/s
%         BALANCE.water.dr_lateralWater = BALANCE.water.dr_lateralWater + (waterflux-excess_water)*1000; % Excess water is removed so that we only keep the net water modification implied by the lateral fluxes
%         fprintf('\t\t\tNet wc change :\t%3.2e m\n',waterflux-excess_water)
%         if excess_water>1e-9
%             GRID.soil.water2pool= GRID.soil.water2pool + excess_water;
%             fprintf('\t\t\tExcess water :\t%3.2e m\n',excess_water)
%             BALANCE.water.dr_lateralExcess=BALANCE.water.dr_lateralExcess + excess_water*1000;            % Added by Leo to have the lateral fluxes in BALANCE
%         end
%         %if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1 || strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoirNoInflow')==1
%         BALANCE.water.dr_DarcyReservoir = BALANCE.water.dr_DarcyReservoir + boundary_water_flux*1000;
%         fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
%         %end
%     end
% 
% end
