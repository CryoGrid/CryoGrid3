% this script is run when the model output is saved (default every year)
% and calculates some "high-level" diagnostics which can then be saved
% separately for each realization and year 
function DIAG = diagnose_output_yearly( OUT, PARA, GRID, FORCING )

    DIAG = {};
    ts = OUT.timestamp();
    dt= abs(ts(1)-ts(2));
    
    %MAAT
    [~, startForcingIndex] = min( abs( ts(1) - FORCING.data.t_span ) );
    [~, endForcingIndex] = min( abs( ts(end) - FORCING.data.t_span ) );
    DIAG.maat = mean( FORCING.data.Tair( startForcingIndex:endForcingIndex ) );
    
    % MAGST, MAGT [°C]
    soil_surface_altitude = OUT.location.soil_altitude();
    altitude_grid = PARA.location.initial_altitude - GRID.general.cT_grid;        
    A = zeros( length(altitude_grid), length(soil_surface_altitude) );
    for j=1:size(A,2)   
        A(:,j) = soil_surface_altitude(j)- altitude_grid;      % distance of each grid cell (first dim) to the surface for each timestep (second dim)
    end
    
    % MAGST as MAGT in 0 m depth
    requestedDepth = 0.01;
    [~, indexes] = min( abs( A - requestedDepth ) );
    linindexes = sub2ind( [length(altitude_grid),length(ts)], indexes, [1:1:length(ts)] );
    DIAG.magst = mean( OUT.cryoGrid3(linindexes) );
    
    % MAGT in 1 m depth
    requestedDepth = 1.0;
    [~, indexes] = min( abs( A - requestedDepth ) );
    linindexes = sub2ind( [length(altitude_grid),length(ts)], indexes, [1:1:length(ts)] );
    DIAG.magt1m = mean( OUT.cryoGrid3(linindexes) );

    % MAGT in 10 m depth
    requestedDepth = 10.0;
    [~, indexes] = min( abs( A - requestedDepth ) );
    linindexes = sub2ind( [length(altitude_grid),length(ts)], indexes, [1:1:length(ts)] );
    DIAG.magt10m = mean( OUT.cryoGrid3(linindexes) );
    
    % MAGT in 30 m depth
    requestedDepth = 30.0;
    [~, indexes] = min( abs( A - requestedDepth ) );
    linindexes = sub2ind( [length(altitude_grid),length(ts)], indexes, [1:1:length(ts)] );
    DIAG.magt30m = mean( OUT.cryoGrid3(linindexes) );    
    
    % maxiumum thaw depth (TDmax) [ m ] % this does not work for very deep
    % thawing, since infiltration altitude is not defined when surface
    % frozen 
    [ ~, idx ] = min( OUT.location.infiltration_altitude() );
    DIAG.TDmax = OUT.location.soil_altitude(idx) - OUT.location.infiltration_altitude(idx) ;
    DIAG.t_TDmax = ts(idx);
    DIAG.doy_TDmax = day( datetime( ts(idx), 'ConvertFrom','datenum'), 'dayofyear' );
    
    % end-of-summer water table depth (WTDeos) [ m ]
    DIAG.WTDeos = OUT.location.water_table_altitude(idx) - OUT.location.soil_altitude(idx);
    
    % end-of-winter snow depth (SDmax) [ m ]
    [ DIAG.SDmax, idx] = nanmax( abs( OUT.location.surface_altitude(1:round(length(ts)/2)) - OUT.location.altitude(1:round(length(ts)/2)) ) );   %search only in first half of year
    DIAG.t_SDmax = ts(idx);
    DIAG.doy_SDmax = day( datetime( ts(idx), 'ConvertFrom','datenum'), 'dayofyear' );

    % start of snow-free season ( = thawing period ) [ UTC Julian day ]
    % note that the thawing period is determined per realization and not
    % for the entire landscape
    startSnowFreeIndex = find ( isnan( OUT.snow.topPosition() ), 1, 'first' );
    DIAG.t_startSnowFree = ceil ( ts( startSnowFreeIndex ) );
    DIAG.doy_startSnowFree = day( datetime( ts(startSnowFreeIndex), 'ConvertFrom','datenum'), 'dayofyear' );

    % end of snow-free season ( = thawing period ) [ UTC Julian day ]
    endSnowFreeIndex = length(ts) + 1 - find ( isnan( flip( OUT.snow.topPosition()) ), 1, 'first' );
    DIAG.t_endSnowFree = floor ( ts( endSnowFreeIndex ) );
    DIAG.doy_endSnowFree = day( datetime( ts(endSnowFreeIndex), 'ConvertFrom','datenum'), 'dayofyear' );

    % mean SEB during thawing period [ W / m^2 ]
    DIAG.meanQnet = sum( OUT.EB.Qnet ( startSnowFreeIndex:endSnowFreeIndex ) .* dt ) ./ abs( ts(startSnowFreeIndex) - ts(endSnowFreeIndex) )  ;
    DIAG.meanQh = sum( OUT.EB.Qh ( startSnowFreeIndex:endSnowFreeIndex ) .* dt ) ./ abs( ts(startSnowFreeIndex) - ts(endSnowFreeIndex) )  ;
    DIAG.meanQe = sum( OUT.EB.Qe ( startSnowFreeIndex:endSnowFreeIndex ) .* dt ) ./ abs( ts(startSnowFreeIndex) - ts(endSnowFreeIndex) )  ;
    DIAG.meanQg = sum( OUT.EB.Qg ( startSnowFreeIndex:endSnowFreeIndex ) .* dt ) ./ abs( ts(startSnowFreeIndex) - ts(endSnowFreeIndex) )  ;
    
    % accumulated water balance during thawing period [ mm ]
    DIAG.accumP = sum( OUT.WB.dp_rain( startSnowFreeIndex:endSnowFreeIndex ) + OUT.WB.dp_snow( startSnowFreeIndex:endSnowFreeIndex ) );
    DIAG.accumET = sum( OUT.WB.de( startSnowFreeIndex:endSnowFreeIndex ) + OUT.WB.ds( startSnowFreeIndex:endSnowFreeIndex ) );
    DIAG.accumM = max( OUT.WB.dW_soil() );  % this is a bit cheated, but should be right
    DIAG.accumRint = sum(   OUT.WB.dr_lateralWater( startSnowFreeIndex:endSnowFreeIndex ) - ...
                            OUT.WB.dr_DarcyReservoir( startSnowFreeIndex:endSnowFreeIndex ) + ...
                            OUT.WB.dr_lateralSnow( startSnowFreeIndex:endSnowFreeIndex ) );      %should be 0
    DIAG.accumRext = sum(   OUT.WB.dr_DarcyReservoir( startSnowFreeIndex:endSnowFreeIndex ) + ...
                            OUT.WB.dr_surface( startSnowFreeIndex:endSnowFreeIndex ) + ...       %should be 0
                            OUT.WB.dr_lateralExcess( startSnowFreeIndex:endSnowFreeIndex ) + ... %should be 0
                            OUT.WB.dr_external( startSnowFreeIndex:endSnowFreeIndex ) + ...      %should be 0
                            OUT.WB.dr_excessSnow( startSnowFreeIndex:endSnowFreeIndex ) );       %should be 0
    DIAG.accumDeltaWsoil = sum( OUT.WB.dW_soil( startSnowFreeIndex:endSnowFreeIndex ) ); 
    DIAG.accumDeltaS = DIAG.accumP + DIAG.accumM + DIAG.accumET + DIAG.accumRint + DIAG.accumRext;
    
    %total excess ice melt [ m ]
    DIAG.excessIceThawed = sum( OUT.xice.excessIceThawed );
    
    %total ground subsidence [ m ]
    DIAG.subsidence = OUT.location.soil_altitude(end) - OUT.location.soil_altitude(1);
    
    %accumulated sediment fluxes per worker (no necessarily applied yet)
    DIAG.accumSedimentOrganic = sum( sum( OUT.lateral.sediment_fluxes_o ) );
    DIAG.accumSedimentMineral = sum( sum( OUT.lateral.sediment_fluxes_m ) );
    DIAG.accumSedimentDiff = sum( sum( OUT.lateral.sediment_fluxes_diff ) );
    DIAG.accumSedimentAdv = sum( sum( OUT.lateral.sediment_fluxes_adv ) );
