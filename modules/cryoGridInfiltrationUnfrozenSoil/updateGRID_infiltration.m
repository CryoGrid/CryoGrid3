function [ wc, GRID, surface_runoff ] = updateGRID_infiltration(wc, GRID, PARA, surface_runoff)

    %%% step 2: GRID update


    soilGRIDsizeOld = sum(GRID.soil.cT_domain);

    %%% step 2a) remove cells filled with air (e.g. due to evaporation
    %%% of uppermost grid cell )
    while (GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)+wc(1)<1e-6)
        disp('infiltration - update GRID - removing air cell') %4

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
            PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub) < (PARA.location.absolute_maxWater_altitude -1e-9)  %tsvd fix for ...
        %PARA.location.initial_altitude-GRID.general.K_grid(GRID.soil.cT_domain_ub)<PARA.location.absolute_maxWater_altitude

        %wc(1)>=1   % this prevents a bug for very small
            %surface_runoff when upper cell not filled // but this
            %does not allow ponding on top of actual soil
        disp('infiltration - update GRID - ponding of water below water table') %1

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
        GRID.soil.cT_actPor =   [ 1; GRID.soil.cT_actPor ];                         % set actual porosity to 1
        GRID.soil.cT_mineral =  [ 0 ; GRID.soil.cT_mineral ];
        GRID.soil.cT_soilType = [ 3; GRID.soil.cT_soilType];                        % soilType 3 = pond (sand freeze curve, field cap =0 )

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
    cellsChanged = soilGRIDsizeNew - soilGRIDsizeOld;
    if cellsChanged > 0
        disp( [ 'infiltration - reinitializing LUT - ', num2str(cellsChanged), ' new water cell(s)' ] ); %2
        GRID.soil.cT_water = wc;
        GRID = initializeSoilThermalProperties(GRID, PARA);   
    elseif cellsChanged < 0
        disp('infiltration - shortening LUT - removed water cell(s)');
        GRID.soil.cT_water(1:-cellsChanged) = [];
        GRID.soil.cT_frozen(1:-cellsChanged) = [];
        GRID.soil.cT_thawed(1:-cellsChanged) = [];
        GRID.soil.K_frozen(1:-cellsChanged) = [];
        GRID.soil.K_thawed(1:-cellsChanged) = [];
        GRID.soil.conductivity(1:-cellsChanged,:) = [];
        GRID.soil.capacity(1:-cellsChanged,:) = [];
        GRID.soil.liquidWaterContent(1:-cellsChanged,:) = [];
    end
    
    % step 2d) update GRID domains of water body / lake
    if GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-6    % upper soil cell pure air/water
        % general water body extent
        cT_waterBody = GRID.soil.cT_organic+GRID.soil.cT_mineral<=1e-6;
        GRID.lake.cT_domain(logical(GRID.air.cT_domain+GRID.snow.cT_domain)) = 0;
        GRID.lake.cT_domain(GRID.soil.cT_domain) = cT_waterBody;
        [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
        GRID.lake.K_domain(logical(GRID.air.K_domain+GRID.snow.K_domain)) = 0;
        GRID.lake.K_domain(GRID.lake.cT_domain_ub:GRID.lake.cT_domain_lb+1) = 1;
        GRID.lake.K_domain(GRID.lake.cT_domain_lb+2:end) = 0;
        [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub] = LayerIndex(GRID.lake.K_domain);
    else
        GRID.lake.cT_domain = false(size(GRID.general.cT_grid));
        GRID.lake.K_domain = false(size(GRID.general.K_grid));
        [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
        [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub] = LayerIndex(GRID.lake.K_domain);
    end



end