function [T, GRID, FLAKE]=cryoGridFLAKE(FLAKE,GRID,SEB,PARA,T,T_old,timestep,k_eff,dE_dt,dE_dt_SEB)

%update average temperature of the water domain
FLAKE.t_mnw_n_flk = sum(T_old(GRID.lake.water.cT_domain).*GRID.general.K_delta(GRID.lake.water.cT_domain)) ...
                 ./ sum(GRID.general.K_delta(GRID.lake.water.cT_domain)) + 273.15;
                                                                    
%update water depth
FLAKE.depth_w = GRID.general.K_grid(GRID.lake.water.cT_domain_lb + 1) - GRID.general.K_grid(GRID.lake.water.cT_domain_ub);
      
%update mean temperature according to changes in average temperature with
%water depth. This prevents changes in the bottom temperature due to sudden
%changes in water depth 
if GRID.lake.ice.z_ice>0
    FLAKE.t_mnw_p_flk = FLAKE.t_mnw_n_flk;
    FLAKE.t_mnw_n_flk = FLAKE.t_wml_n_flk - FLAKE.c_t_n_flk.*(FLAKE.t_wml_n_flk-FLAKE.t_bot_n_flk).*(1.-FLAKE.h_ml_n_flk./FLAKE.depth_w);    
    GRID.lake.ice.z_ice = GRID.lake.ice.z_ice - (FLAKE.t_mnw_p_flk-FLAKE.t_mnw_n_flk).*4.2e6.*FLAKE.depth_w./(334e3 * 910);
end

%update ice cover status and thickness
if FLAKE.h_ice_n_flk<1e-3 && GRID.lake.ice.z_ice>=1e-3
    FLAKE.l_ice_create = true;
else
    FLAKE.l_ice_create = false;
end
FLAKE.h_ice_n_flk = GRID.lake.ice.z_ice;

%set boundary conditions
FLAKE.q_snow_flk= 0.;
FLAKE.q_ice_flk = 0.;
FLAKE.q_bot_flk = -(T_old(GRID.lake.water.cT_domain_lb+1)-T_old(GRID.lake.water.cT_domain_lb))...
               ./ GRID.general.cT_delta(GRID.lake.water.cT_domain_lb) .* k_eff(GRID.lake.water.cT_domain_lb+1);

FLAKE.q_w_flk = sum(dE_dt_SEB(GRID.lake.water.cT_domain_ub:GRID.lake.water.cT_domain_lb+1)) ...
              - SEB.Sin_water ...
              + sum(dE_dt(GRID.lake.water.cT_domain_ub:GRID.lake.water.cT_domain_lb)) ...
              - sum(dE_dt_SEB(GRID.lake.water.cT_domain_ub:GRID.lake.water.cT_domain_lb)) ...
              + FLAKE.q_bot_flk ...
              + GRID.lake.ice.dz_dt_freeze*334e3*910 ...
              + GRID.lake.ice.dE_dt_melt_residual;

if FLAKE.q_w_flk>0 && GRID.lake.ice.z_ice>0
    %correct for positive water heat fluxes under ice cover. This is
    %necessary due to the very simple ice cover scheme which needs to
    %be improved in fututre model versions.    
    GRID.lake.ice.z_ice = GRID.lake.ice.z_ice - FLAKE.q_w_flk.*(timestep *24*3600)/(334e3 * 910);
    FLAKE.q_w_flk = 0;
    FLAKE.h_ice_n_flk = GRID.lake.ice.z_ice;
end

%transform u_star in atmosphere to u_star in water assuming constant
%desities for air and water.
FLAKE.u_star_w_flk=(SEB.u_star.^2.*(1.293./1e3)).^0.5;

%calculate FLAKE radiative heat transfer
[FLAKE.i_atm_flk,...
 FLAKE.i_w_flk,...
 FLAKE.i_ice_flk,...
 FLAKE.i_snow_flk,...
 FLAKE.i_h_flk,...
 FLAKE.i_bot_flk,...
 FLAKE.i_intm_0_h_flk,...
 FLAKE.i_intm_h_d_flk] = flake_radflux(SEB.Sin_water,FLAKE,PARA);

% %call FLAKE main function
% if PARA.technical.MEX
%     FLAKE = flake_driver_mex(FLAKE, timestep.*24.*3600);
% else
    FLAKE = flake_driver(FLAKE, timestep.*24.*3600);
% end

%map water temperature from shape function on the regular grid
%calculate dimensionless water depth of thermocline
% zeta_cT = (GRID.lake.water.cT_grid(GRID.lake.water.K_grid>=FLAKE.h_ml_n_flk) - FLAKE.h_ml_n_flk) ...
%        ./ ((FLAKE.depth_w - FLAKE.h_ml_n_flk));

zeta_K  = (GRID.lake.water.K_grid(GRID.lake.water.K_grid>=FLAKE.h_ml_n_flk) - FLAKE.h_ml_n_flk) ...
       ./ ((FLAKE.depth_w - FLAKE.h_ml_n_flk));
zeta_K  = [zeta_K; 1];
zeta_delta = diff(zeta_K);

if ~isempty(zeta_delta)
    %The original FLAKE shape function not used due to interpolation
    %errors on discrete grid
    %phi_theta = (40/3*FLAKE.c_t_n_flk-20/3)*zeta_cT ...
    %          + (18-30*FLAKE.c_t_n_flk)*zeta_cT.^2 ...
    %          + (20*FLAKE.c_t_n_flk-12)*zeta_cT.^3 ...
    %          + (5/3 - 10/3*FLAKE.c_t_n_flk)*zeta_cT.^4;
    
    %Integrated shape function to calculate averages on discrete grid
    int_phi_theta = (40/3*FLAKE.c_t_n_flk-20/3)*1/2*zeta_K.^2 ...
                  + (18-30*FLAKE.c_t_n_flk)*1/3*zeta_K.^3 ...
                  + (20*FLAKE.c_t_n_flk-12)*1/4*zeta_K.^4 ...
                  + (5/3 - 10/3*FLAKE.c_t_n_flk)*1/5*zeta_K.^5;
    
    %calculate partial intergral for grid cells
    int_phi_theta = int_phi_theta(2:end) - int_phi_theta(1:end-1);
    phi_theta     = int_phi_theta./zeta_delta;
    
    %change from dimensionless temperature to absolute temperature
    Tw_z = -(phi_theta*(FLAKE.t_wml_n_flk-FLAKE.t_bot_n_flk) - FLAKE.t_wml_n_flk)-273.15;
else
    Tw_z = [];
end
%map temperature on water grid
Tw = T(GRID.lake.water.cT_domain);

%note that a step is introduced which introduces small discretization
%errors which need to be corrected later
Tw(GRID.lake.water.K_grid>=FLAKE.h_ml_n_flk) = Tw_z;
Tw(GRID.lake.water.K_grid<FLAKE.h_ml_n_flk)  = FLAKE.t_wml_n_flk-273.15;
T(GRID.lake.water.cT_domain) = Tw;

%calculate difference between average FLAKE temperature and gridded
%temperature. Minor interpolation errors occure due to discretesation
%of the mixed layer depth. In adiition, errors can occur due to FLAKE
%constrains which do not account for convection under ice due to
%e.g. radiative heating. The missing energy is assumed to finally
%accumulate in ice melt. Thus, temperature difference is corrected
%with ice thickness.
dT = mean(Tw)-((FLAKE.t_mnw_n_flk-273.15));
dE = dT.*4.2e6.*FLAKE.depth_w;
if GRID.lake.ice.z_ice>0
    GRID.lake.ice.z_ice = GRID.lake.ice.z_ice + dE./(334e3 * 910);   
    %set T water surf to zero
    %dE=(0-T(GRID.lake.water.cT_domain_ub))*GRID.general.K_delta(GRID.lake.water.cT_domain_ub)*4.2e6;
    %GRID.lake.ice.z_ice = GRID.lake.ice.z_ice + dE./(334e3 * 910);
    %T(GRID.lake.water.cT_domain_ub) = 0;
else
    ub_h_ml=min([GRID.lake.water.K_grid(GRID.lake.water.K_grid>=FLAKE.h_ml_n_flk); FLAKE.depth_w]);
    lb_h_ml=FLAKE.depth_w;
    if (lb_h_ml-ub_h_ml)>0
        Tw(GRID.lake.water.K_grid>=FLAKE.h_ml_n_flk) = Tw_z - dE/(4.2e6* (lb_h_ml-ub_h_ml));
    else
        Tw = Tw - dE/(4.2e6 * FLAKE.depth_w);
    end
    T(GRID.lake.water.cT_domain) = Tw;
end

%update FLAKE average temperature of the water domain
FLAKE.t_mnw_n_flk = sum(T(GRID.lake.water.cT_domain).*GRID.general.K_delta(GRID.lake.water.cT_domain)) ...
                 ./ sum(GRID.general.K_delta(GRID.lake.water.cT_domain)) + 273.15;
                         
             
             
