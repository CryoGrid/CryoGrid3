% Plotting script for temperature and water content fields
% Author: Jan Nitzbon
% 
%function plot_output(dirname, runname, number_of_realizations)

%clear all
close all

 dirname = '../runs/';
 runname = 'SPINUP-EXICE_197906-201406_stratSamExice_rf1_sf1_maxSnow0.40_snowDens=200_wt0.0_extFlux0.0000_fc0.30_exice0.90_natPor0.40';
 year = 1982;
 number_of_realizations = 1;
% 

    %% load output data and settings from files of all workers


    dir = dirname;%'/home/jnitzbon/gls1/CryoGrid/CryoGrid3_infiltration_xice_mpi_DEV/';
    cm=load( [ './cm_blueautumn.mat' ] );

    infil=0;
    xice=0;
    rainFrac=0;
    snowFrac=0;
    wt = 0.0;

    OUTS = {};

    for i=1:number_of_realizations
        run = runname; %sprintf('testrun_mpi_waterExchange_noGroundHeatFlux_realization' );
        if number_of_realizations>1
            run = [ run num2str(i) ];
        end
        outputfile = [dir run  '/' run '_output' num2str(year) '.mat'];
        configfile = [dir run  '/' run '_settings.mat'];

        load(outputfile);
        load(configfile);

        OUTS{i}.OUT=OUT;
        OUTS{i}.PARA = PARA;
        OUTS{i}.GRID = GRID;
        OUTS{i}.FORCING = FORCING;

        clear OUT PARA GRID FORCING
    end
    
    %% extract relevant vectors

    ts = OUTS{1}.OUT.timestamp();
    zs = {};
    Ts = {};
    LWCs = {};
    WCs = {};
    snowTop = {};
    soilTop = {};

    for i=1:number_of_realizations
        %zs{i} = OUTS{i}.PARA.ensemble.initial_altitude(i)-OUTS{i}.GRID.general.cT_grid;
        zs{i} = OUTS{i}.PARA.location.altitude-OUTS{i}.GRID.general.cT_grid;        
        Ts{i} = OUTS{i}.OUT.cryoGrid3;
        LWCs{i} = OUTS{i}.OUT.liquidWater;
        WCs{i} = OUTS{i}.OUT.water;
        snowTop{i} = OUTS{i}.OUT.soil.topPosition;
    end
    %lakeFloor = OUT.soil.lakeFloor();
    %lakeFloor = [ NaN(length(soilTop)-length(lakeFloor),1); lakeFloor ];

    % limits
    minz = min(OUTS{1}.PARA.location.altitude - 10);
    maxz = max(OUTS{1}.PARA.location.altitude + 0.5);

    mint = min(ts);
    maxt = max(ts);


    %% plot evolution of simple ensemlbe variables

    if number_of_realizations>1
        
        figure;

        subplot(2,3,1)
        hold on
        for i=1:number_of_realizations
            plot(ts, OUTS{1}.OUT.ensemble.active_layer_depth_altitude(:,i)) ; %abs(OUTS{1}.OUT.ensemble.altitude(:,i)-
        end
        datetick;
        title('Frost table [m asl]')
        hold off

        subplot(2,3,2)
        hold on
        for i=1:number_of_realizations
            plot(ts, OUTS{1}.OUT.ensemble.water_table(:,i));
        end
        datetick;
        title('Water table [m asl]')
        hold off

        subplot(2,3,3)
        hold on
        for i=1:number_of_realizations
            plot(ts, OUTS{1}.OUT.ensemble.surface_altitude(:,i));
        end
        datetick;
        title('Surface altitude [m asl]')
        hold off

        subplot(2,3,4)
        hold on
        for i=1:number_of_realizations
            plot(ts, nansum(OUTS{i}.OUT.ensemble.heat_fluxes , 2));
        end
        datetick;
        title('Heat fluxes [ W / m^2 ]')
        hold off

        subplot(2,3,5)
        hold on
        for i=1:number_of_realizations
            plot(ts, nansum(OUTS{i}.OUT.ensemble.water_fluxes .*1000 .* 3600 , 2));
        end
        datetick;
        title('Water fluxes [mm/h]')
        hold off

        subplot(2,3,6)
        hold on
        for i=1:number_of_realizations
            plot(ts, nansum(OUTS{i}.OUT.ensemble.snow_fluxes .*1000 .* 3600 , 2) );
        end
        datetick;
        title('Snow fluxes [mm/h]')
        hold off

        currentFigure = gcf;
        title(currentFigure.Children(end), runname);
        
    end

        

    %% plot temperature and water content fields

    figure;

    title(runname);

    %temperature fields
    for i=1:number_of_realizations
        subplot(2, number_of_realizations, i);
        pcolor( ts', zs{i}', Ts{i});
        hold on;
        shading flat; % do not show grid
        %shading interp;
        % colormap and colorbar
        caxis( [ -40, 20] );
        colormap(gca, cm.Colormap_blueautumn);
        if i==1
            cbar=colorbar('location','westoutside');
        end
        % layout
        axis( [ mint maxt minz maxz ] );
        datetick;
        xlabel('time')
        ylabel('$z$ [m]', 'Interpreter', 'latex');
        xlabel(cbar, '$T$ [$^\circ$C]', 'Interpreter', 'latex');
        title(sprintf('Temperature - realization %d', i));
        hold off;
    end
    
    % (liquid) water content fields
    for i=1:number_of_realizations
        subplot(2, number_of_realizations, number_of_realizations+i);
        pcolor( ts', zs{i}', WCs{i});  %ax,
        hold on;
        caxis( [ 0. , 1. ] );
        colormap(gca, 'parula');
        colormap(gca, flipud(colormap));
        if i==1
            cbar=colorbar('location','westoutside');
        end
        shading flat;
        % layout
        axis( [ mint maxt minz maxz ] );
        datetick;
        xlabel('time');
        ylabel('$z$ [m]', 'Interpreter', 'latex');
        xlabel(cbar, '$\theta_w$ [-]', 'Interpreter', 'latex');
        title(sprintf('Water content - realization %d', i));
        hold off;
    end
    
    currentFigure = gcf;
    title(currentFigure.Children(end), runname);    
    
    
    %% plot soil moisture and temperature of requested depth
    
    % request time series for a certain depths
    requestedDepth =  0.3 ;
    
    % get the altitude of the uppermost soil cell (which indeed contains
    % mineral or organic material, i.e. excluding a potential waterbody)
    soil_surface_altitude = OUTS{1}.PARA.location.altitude + min( OUTS{1}.OUT.soil.topPosition, OUTS{1}.OUT.soil.lakeFloor );
    
    % the static altitude grid
    altitude_grid = zs{1};
    
    % compute the index of the requested cell for each timestep
    A = zeros( length(altitude_grid), length(soil_surface_altitude) );   % matrix to serach in
    for j=1:size(A,2)   %loop over all timesteps       
        A(:,j) = soil_surface_altitude(j)- altitude_grid;      % distance of each grid cell (first dim) to the surface for each timestep (second dim)
    end   
    [~, indexes] = min( abs( A - requestedDepth ) );    % determine index of closest cell to the requested depth
    
    % transform column-wise index to linear index of whole matrix
    linindexes = sub2ind( [length(altitude_grid),length(ts)], indexes, [1:1:length(ts)] );

    figure;
    plot( ts, LWCs{i}(linindexes) );

    figure;
    plot( ts, Ts{i}(linindexes) );