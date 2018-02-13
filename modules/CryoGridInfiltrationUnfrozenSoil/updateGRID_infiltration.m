function [ wc, GRID, surface_runoff ] = updateGRID_infiltration(wc, GRID, PARA, surface_runoff)

    %%% step 2: GRID update
    %%% TODO: add a function updateGRID_infiltration

    soilGRIDsizeOld = sum(GRID.soil.cT_domain);

    %%% step 2a) remove cells filled with air (e.g. due to evaporation
    %%% of uppermost grid cell )
    while (GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)+wc(1)<=0)
        disp('infiltration - update GRID - removing air cell')

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

        GRID.soil.cT_organic(1)=[];
        GRID.soil.cT_natPor(1)=[];
        GRID.soil.cT_actPor(1)=[];
        GRID.soil.cT_mineral(1)=[];
        GRID.soil.cT_soilType(1)=[];


        GRID.soil.excessGroundIce(1)=[];

    end

    %%% step 2b) ponding of surface runoff below water table
    while surface_runoff>1e-6 && ...                                % not >0 as sometimes numerical errors occur during calculation of surface_runoff
            PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub)<PARA.location.absolute_maxWater_altitude
            %wc(1)>=1   % this prevents a bug for very small
            %surface_runoff when upper cell not filled // but this
            %does not allow ponding on top of actual soil
        disp('infiltration - update GRID - ponding of water below water table')

        h = PARA.location.absolute_maxWater_altitude - ( PARA.location.initial_altitude - GRID.general.K_grid(GRID.soil.cT_domain_ub) ) ;   % this is guruanteed to be >0

        % create new water cell / change GRID domains
        GRID.soil.cT_domain(GRID.air.cT_domain_lb)=1;
        GRID.soil.K_domain(GRID.air.K_domain_lb)=1;
        GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub-1;
        GRID.soil.K_domain_ub=GRID.soil.K_domain_ub-1;
        GRID.air.cT_domain(GRID.air.cT_domain_lb)=0;
        GRID.air.K_domain(GRID.air.K_domain_lb)=0;
        GRID.air.cT_domain_lb=GRID.air.cT_domain_lb-1;
        GRID.air.K_domain_lb=GRID.air.K_domain_lb-1;

        % fill new water cell
        cellSize = PARA.technical.waterCellSize;
        waterAdded = min( [surface_runoff, cellSize, h] ); % add water until water table is reached or surface_runoff "empty"
        wc = [ waterAdded./cellSize ; wc ];
        surface_runoff = surface_runoff - waterAdded;

        % update remaining soil fields with exception of cT_water
        GRID.soil.cT_organic =  [ 0 ; GRID.soil.cT_organic ];
        GRID.soil.cT_natPor =   [ GRID.soil.cT_natPor(1); GRID.soil.cT_natPor ];    % take natPor of cell below
        GRID.soil.cT_actPor =   [ 1; GRID.soil.cT_actPor ];                         % set to 1
        GRID.soil.cT_mineral =  [ 0 ; GRID.soil.cT_mineral ];
        GRID.soil.cT_soilType = [ 1; GRID.soil.cT_soilType];                        % assume sand as soil type for water cell

        GRID.soil.excessGroundIce = [ 0 ; GRID.soil.excessGroundIce ];

        % update GRID spacings
        GRID.general.K_grid(GRID.soil.cT_domain_ub) = GRID.general.K_grid(GRID.soil.cT_domain_ub+1)-cellSize;
        GRID.general.K_grid(GRID.air.cT_domain) = [GRID.general.K_grid(GRID.air.cT_domain_lb)+(-GRID.snow.snowCellSize)*(GRID.air.cT_domain_lb-1):GRID.snow.snowCellSize:GRID.general.K_grid(GRID.air.cT_domain_lb)]';
        GRID.general.cT_grid = ( GRID.general.K_grid(1:end-1)+ GRID.general.K_grid(2:end))/2; %grid on which capacity and temperature information lives (midpoints of grid cells)
        GRID.general.cT_delta = (- GRID.general.cT_grid(1:end-1,1)+ GRID.general.cT_grid(2:end,1));
        GRID.general.K_delta = (- GRID.general.K_grid(1:end-1,1)+ GRID.general.K_grid(2:end,1));
        GRID.soil.soilGrid = [ GRID.general.K_grid(GRID.soil.cT_domain_ub) ; GRID.soil.soilGrid ];

    end
    
    
    
    %%% step 2c)  check if soil/air domains changed --> LUT update
    soilGRIDsizeNew = sum(GRID.soil.cT_domain);
    if soilGRIDsizeOld~=soilGRIDsizeNew
        disp('infiltration - reinitializing LUT - soil/air domains changed');
        GRID.soil.cT_water = wc;
        GRID = initializeSoilThermalProperties(GRID, PARA);   
    end



end
