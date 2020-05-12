function [wc, GRID, BALANCE] = CryoGridLateralWater( PARA, GRID, BALANCE, T, wc)

    labBarrier();
    % check preconditions
    precondition_waterExchange = checkPreconditionWaterExchange( T, GRID );
    if precondition_waterExchange
        % WRAPPER
        fprintf('\t\t\tsync - exchanging water\n');
        % calculate lateral water fluxes
        water_fluxes= zeros(numlabs,numlabs); % in m of height change
        PACKAGE_waterExchange.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
        PACKAGE_waterExchange.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
        PACKAGE_waterExchange.infiltration_condition = T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub);
        PACKAGE_waterExchange.soil_altitude = PARA.ensemble.soil_altitude(labindex);

        for j=1:numlabs
            if PARA.ensemble.adjacency_water(labindex,j)==1 % only send between connected workers to save computation time
                labSend( PACKAGE_waterExchange, j, 20);
            end
        end
        for j=1:numlabs
            if PARA.ensemble.adjacency_water(labindex,j)==1
                PACKAGE_waterExchange_j = labReceive(j, 20);
                water_fluxes = calculateLateralSurfaceSubsurfaceWaterDarcyFluxes( T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes);  % matrix containing all fluxes in [m/s] scaled to row index
            end
        end
                
        % Calculate possible boundary fluxes
        [ boundary_water_flux ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T); % this could also be separated for surf/subs contributions
        
        % Check for water availability and set real water fluxes
        [ water_fluxes_worker, boundary_water_flux ] = calculateLateralWaterAvailable( PARA,GRID, wc, water_fluxes, boundary_water_flux );
        water_fluxes_gather=zeros(numlabs,numlabs,numlabs);
        water_fluxes_gather(:,:,labindex)=water_fluxes_worker;

        % Send real fluxes all around
        for j=1:numlabs
            if j~=labindex
                labSend( water_fluxes_worker, j, 22);
            end
        end

        for j=1:numlabs
            if j~=labindex
                water_fluxes_worker_j = labReceive(j, 22);
                water_fluxes_gather(:,:,j)=water_fluxes_worker_j;
            end
        end
        water_fluxes=nansum(water_fluxes_gather,3);
        
        waterflux=nansum(water_fluxes(labindex,:))+boundary_water_flux;

        % apply lateral water flux directly (as bulk subsurface flux)
        [wc, excess_water, lacking_water] = bucketScheme(T, wc, zeros( size(wc) ), GRID, PARA, waterflux);
        try
            assert( lacking_water < 1e-9, 'CryoGrid3 - lateral exchange - lacking water>0');    % there should be no lacking water as this was checked for
        catch
            fprintf('\t\t\tLacking water = %3.2e m\n', lacking_water );
        end

        % Store and display
        BALANCE.water.lateral_water_fluxes_vector=BALANCE.water.lateral_water_fluxes_vector + water_fluxes(labindex,:); % here we decide in which units we want it. now in [m] height change
        BALANCE.water.lateral_water_fluxes_matrix=BALANCE.water.lateral_water_fluxes_matrix + water_fluxes(:,:); % here we decide in which units we want it. now in [m] height change
        BALANCE.water.lateral_water_fluxes_boundary=BALANCE.water.lateral_water_fluxes_boundary + boundary_water_flux;
        
        BALANCE.water.dr_lateralWater = BALANCE.water.dr_lateralWater + (waterflux-excess_water)*1000; % Excess water is removed so that we only keep the net water modification implied by the lateral fluxes
        fprintf('\t\t\tNet wc change :\t%3.2e m\n',waterflux-excess_water)
        if excess_water>1e-9
            GRID.soil.water2pool= GRID.soil.water2pool + excess_water;
            fprintf('\t\t\tExcess water :\t%3.2e m\n',excess_water)
            BALANCE.water.dr_lateralExcess=BALANCE.water.dr_lateralExcess + excess_water*1000;
        end
        BALANCE.water.dr_DarcyReservoir = BALANCE.water.dr_DarcyReservoir + boundary_water_flux*1000;
        fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
    end

end