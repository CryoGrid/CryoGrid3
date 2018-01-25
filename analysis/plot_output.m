% Plotting script for temperature and water content fields
% Author: Jan Nitzbon
% 
function plot_output(dirname, runname, number_of_realizations)

% dirname = './';
% runname = 'testrunMPI_PAIRWISE_xH1_xW1_xS1_infil1_xice1_rF1.0_sF1.0_Dsnow1e-6_realization';
% number_of_realizations = 2;
% 

    %% load output data and settings from files of all workers


    dir = dirname;%'/home/jnitzbon/gls1/CryoGrid/CryoGrid3_infiltration_xice_mpi_DEV/';
    cm=load( [ dir 'cm_blueautumn.mat' ] );

    infil=0;
    xice=0;
    rainFrac=0;
    snowFrac=0;
    wt = 0.0;

    OUTS = {};

    for i=1:number_of_realizations
        run = runname; %sprintf('testrun_mpi_waterExchange_noGroundHeatFlux_realization' );
        run = [ run num2str(i) ];
        outputfile = [dir run  '/' run '_output.mat'];
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
        zs{i} = OUTS{i}.PARA.ensemble.initial_altitude(i)-OUTS{i}.GRID.general.cT_grid;
        Ts{i} = OUTS{i}.OUT.cryoGrid3;
        LWCs{i} = OUTS{i}.OUT.liquidWater;
        WCs{i} = OUTS{i}.OUT.water;
        snowTop{i} = OUTS{i}.OUT.soil.topPosition;
    end
    %lakeFloor = OUT.soil.lakeFloor();
    %lakeFloor = [ NaN(length(soilTop)-length(lakeFloor),1); lakeFloor ];

    % limits
    minz = min(OUTS{1}.PARA.ensemble.initial_altitude - 4);
    maxz = max(OUTS{1}.PARA.ensemble.initial_altitude + 0.5);

    mint = min(ts);
    maxt = max(ts);


    %% plot evolution of simple ensemlbe variables

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
        pcolor( ts', zs{i}', LWCs{i});  %ax,
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

    %% plot temperature in 2m depth
    
    figure;
    hold on;
    
    title(runname);
    
    depth = 2.0;
    [val, idx] = min( abs(depth-OUTS{1}.GRID.general.cT_grid) );
    
    plot( ts', Ts{1}( idx , :) );
    plot( ts', Ts{2}( idx , :) );
    
    axis( [mint maxt -20 0 ] );
    
    title(sprintf('Temperature in [°C] at z=%f', OUTS{1}.GRID.general.cT_grid(idx) ));
    
    hold off;
    
    
    