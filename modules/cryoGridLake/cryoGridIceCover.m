function [T,GRID,FLAKE]=cryoGridIceCover(GRID,SEB,FLAKE,T,c_temp,timestep)

%construct water grids relative to the water surface
if ~isempty(GRID.lake.water.cT_domain_ub)
    GRID.lake.water.cT_grid = GRID.general.cT_grid(GRID.lake.water.cT_domain) - GRID.general.K_grid(GRID.lake.water.cT_domain_ub);
    GRID.lake.water.K_grid =  GRID.general.K_grid(GRID.lake.water.K_domain)   - GRID.general.K_grid(GRID.lake.water.cT_domain_ub);
else
    GRID.lake.water.cT_grid = [];
    GRID.lake.water.K_grid = [];
end

%initialisize
GRID.lake.ice.dz_dt_freeze=0;
GRID.lake.ice.dz_dt_melt=0;
GRID.lake.ice.dE_dt_melt_residual=0;

%ice cover build up before
if min(T(GRID.lake.water.cT_domain))<0
    
    %distribute heat into the mixed water layer in order to avoid unrealistic freezing of the first water cell
    if FLAKE.h_ml_n_flk>0
        Tw=T(GRID.lake.water.cT_domain);
        K_delta=GRID.general.K_delta(GRID.lake.water.cT_domain);
        
        Tw(GRID.lake.water.K_grid<FLAKE.h_ml_n_flk) = sum(Tw(GRID.lake.water.K_grid<FLAKE.h_ml_n_flk)...
                                              .* K_delta(GRID.lake.water.K_grid<FLAKE.h_ml_n_flk))...
                                              ./ sum(K_delta(GRID.lake.water.K_grid<FLAKE.h_ml_n_flk));
        
        T(GRID.lake.water.cT_domain)=Tw;
    end
    %calculate energy used for ice cover formation 
    dE_dt_ice  = (0-T(GRID.lake.water.cT_domain)) ...
              .* c_temp(GRID.lake.water.cT_domain) ...
              .* GRID.general.K_delta(GRID.lake.water.cT_domain) ...
              ./ (timestep.*24.*3600);
    
    dE_dt_ice = dE_dt_ice .* (T(GRID.lake.water.cT_domain)<0);
    dE_dt_ice = sum(dE_dt_ice);
    
    %set temperature T to 0°C for freezing cells
    T(GRID.lake.water.cT_domain & T<0)=0;
    GRID.lake.ice.dz_dt_freeze = dE_dt_ice ./ (334e3 * 910);
    GRID.lake.ice.dz_dt_freeze(isinf(GRID.lake.ice.dz_dt_freeze) | isnan(GRID.lake.ice.dz_dt_freeze))=0;    
end

%ice cover melt
%set ice cover melt flag to false
GRID.lake.ice.melt_flag=false;
if ~isempty(GRID.lake.ice.cT_domain_ub) && max(T(GRID.lake.ice.cT_domain))>0     
    dE_dt_ice = (0-T(GRID.lake.ice.cT_domain)) ...
             .* c_temp(GRID.lake.ice.cT_domain) ...
             .* GRID.general.K_delta(GRID.lake.ice.cT_domain) ...
             ./ (timestep.*24.*3600);
    
    dE_dt_ice = dE_dt_ice .* (T(GRID.lake.ice.cT_domain)>0);
    dE_dt_ice = sum(dE_dt_ice);
    
    %set temperature T to 0°C for cells with melting ice
    T((GRID.lake.ice.cT_domain) & T>0) = 0;
    GRID.lake.ice.dz_dt_melt = dE_dt_ice./ (334e3 * 910);
    GRID.lake.ice.dz_dt_melt(isinf(GRID.lake.ice.dz_dt_melt) | isnan(GRID.lake.ice.dz_dt_melt)) = 0;
    
    %set melt flag true if the entire ice cover is affected by melting
    if sum(T(GRID.lake.ice.cT_domain))==0;       
        GRID.lake.ice.melt_flag=true;
    end
end

%melt remaining ice cover
GRID.lake.ice.dE_dt_melt_residual = 0;
if isempty(GRID.lake.ice.cT_domain_ub) && GRID.lake.ice.z_ice>0 && T(GRID.lake.water.cT_domain_ub)>0
    GRID.lake.ice.dz_dt_melt = max(-GRID.lake.ice.z_ice/(timestep*24*3600), -SEB.Qg/(334e3 * 910));
    GRID.lake.ice.dE_dt_melt_residual = GRID.lake.ice.dz_dt_melt*(334e3 * 910);
end

%calculate change of ice cover depth
GRID.lake.ice.dz_dt_ice = GRID.lake.ice.dz_dt_freeze + GRID.lake.ice.dz_dt_melt - GRID.lake.ice.dE_dt_cond_residual./(334e3 * 910);

%change ice cover thickness
GRID.lake.ice.z_ice = GRID.lake.ice.z_ice + GRID.lake.ice.dz_dt_ice.*(timestep.*24.*3600);
GRID.lake.ice.z_ice = GRID.lake.ice.z_ice*(GRID.lake.ice.z_ice>0);

%sort water cells under ice cover according to density
if ~isempty(GRID.lake.ice.cT_domain_ub) && ~isempty(GRID.lake.water.cT_domain_ub)
    DW = waterDensity(T(GRID.lake.water.cT_domain));
    [~, XI]=sort(DW);
    DW=T(GRID.lake.water.cT_domain);
    T(GRID.lake.water.cT_domain)=DW(XI);
end






