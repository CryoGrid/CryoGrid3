%------------------------------------------------------------------------------

function [FLAKE]=flake_driver(FLAKE, del_time) %#codegen

%------------------------------------------------------------------------------
%
% Description:
%
%  The main driving routine of the lake model FLake
%  where computations are performed.
%  Advances the surface temperature
%  and other FLake variables one time step.
%  At the moment, the Euler explicit scheme is used.
%
%  Lines embraced with '!_tmp' contain temporary parts of the code.
%  Lines embraced/marked with '!_dev' may be replaced
%  as improved parameterizations are developed and tested.
%  Lines embraced/marked with '!_dm' are DM's comments
%  that may be helpful to a user.
%  Lines embraced/marked with '!_dbg' are used
%  for debugging purposes only.
%
%
% Current Code Owner: DWD, Dmitrii Mironov
%  Phone:  +49-69-8062 2705
%  Fax:    +49-69-8062 3721
%  E-mail: dmitrii.mironov@dwd.de
%------------------------------------------------------------------------------
%
%  Dimensionless constants 
%  in the equations for the mixed-layer depth 
%  and for the shape factor with respect to the temperature profile in the thermocline
  c_cbl_1       = 0.17            ; % Constant in the CBL entrainment equation
  c_cbl_2       = 1.              ; % Constant in the CBL entrainment equation
  c_sbl_zm_n    = 0.5             ; % Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_zm_s    = 10.             ; % Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_zm_i    = 20.             ; % Constant in the ZM1996 equation for the equilibrium SBL depth
  c_relax_h     = 0.030           ; % Constant in the relaxation equation for the SBL depth
  c_relax_c     = 0.0030          ; % Constant in the relaxation equation for the shape factor
                                    % with respect to the temperature profile in the thermocline

%  Parameters of the shape functions 
%  Indices refer to T - thermocline, S - snow, I - ice,
%  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
%  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
%  at "zeta=0" ad "zeta=1", respectively.
  c_t_min       = 0.5      ; % Minimum value of the shape factor C_T (thermocline)
  c_t_max       = 0.8      ; % Maximum value of the shape factor C_T (thermocline)
%  phi_t_pr0_1   = 40./3.   ; % Constant in the expression for the T shape-function derivative 
%  phi_t_pr0_2   = 20./3.   ; % Constant in the expression for the T shape-function derivative 
  c_tt_1        = 11./18.  ; % Constant in the expression for C_TT (thermocline)
  c_tt_2        = 7./45.   ; % Constant in the expression for C_TT (thermocline)
%  c_b1          = 2./3.    ; % Shape factor (upper layer of bottom sediments)
%  c_b2          = 3./5.    ; % Shape factor (lower layer of bottom sediments)
%  phi_b1_pr0    = 2.       ; % B1 shape-function derivative 
%  c_s_lin       = 0.5      ; % Shape factor (linear temperature profile in the snow layer)
%  phi_s_pr0_lin = 1.       ; % S shape-function derivative (linear profile) 
%  c_i_lin       = 0.5      ; % Shape factor (linear temperature profile in the ice layer)
%  phi_i_pr0_lin = 1.       ; % I shape-function derivative (linear profile) 
%  phi_i_pr1_lin = 1.       ; % I shape-function derivative (linear profile) 
%  phi_i_ast_mr  = 2.       ; % Constant in the MR2004 expression for I shape factor
%  c_i_mr        = 1./12.   ; % Constant in the MR2004 expression for I shape factor
%  h_ice_max     = 3.       ; % Maximum ice tickness in 
                             % the Mironov and Ritter (2004, MR2004) ice model [m] 

% Security constants
  h_snow_min_flk = 1.0E-5         ; % Minimum snow thickness [m]
  h_ice_min_flk  = 1.0E-3         ; % Minimum ice thickness [m]
  h_ml_min_flk   = 1.0E-2         ; % Minimum mixed-layer depth [m]
  h_ml_max_flk   = 1.0E+3         ; % Maximum mixed-layer depth [m]
%  h_b1_min_flk   = 1.0E-3         ; % Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6         ; % Minimum value of the surface friction velocity [m s^{-1}]

% Security constant(s)
  c_small_flk    = 1.0E-10        ; % A small number

% Thermodynamic parameters
  tpl_grav          = 9.81       ; % Acceleration due to gravity [m s^{-2}]
  tpl_t_r           = 277.13     ; % Temperature of maximum density of fresh water [K]
  tpl_t_f           = 273.15     ; % Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05 ; % Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03    ; % Maximum density of fresh water [kg m^{-3}]
%  tpl_rho_i         = 9.1E+02    ; % Density of ice [kg m^{-3}]
%  tpl_rho_s_min     = 1.0E+02    ; % Minimum snow density [kg m^{-3}]
%  tpl_rho_s_max     = 4.0E+02    ; % Maximum snow density [kg m^{-3}]
%  tpl_gamma_rho_s   = 2.0E+02    ; % Empirical parameter [kg m^{-4}]  
                                   % in the expression for the snow density 
%  tpl_l_f           = 3.34E+05   ; % Latent heat of fusion [J kg^{-1}]
  tpl_c_w           = 4.2E+03    ; % Specific heat of water [J kg^{-1} K^{-1}]
%  tpl_c_i           = 1.9E+03    ; % Specific heat of ice [J kg^{-1} K^{-1}]
%  tpl_c_s           = 2.1E+03    ; % Specific heat of snow [J kg^{-1} K^{-1}]
%  tpl_kappa_w       = 0.546      ; % Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
%  tpl_kappa_i       = 2.29       ; % Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]


%==============================================================================

%==============================================================================
%
% Declarations

% Input (procedure arguments)

% The lake depth [m]
% real( ireals), intent(in) :: ;
% Depth of the thermally active layer of bottom sediments [m]
%depth_bs                          ,;
% Temperature at the outer edge of
%t_bs                              ,;
% the thermally active layer of bottom sediments [K]
% The Coriolis parameter [s^{-1}]
%par_coriolis;                     
% 'Typical' extinction coefficient of the lake water [m^{-1}],
%extincoef_water_typ               ,;
% used to compute the equilibrium CBL depth
% The model time step [s]
%del_time                          ,;
% Surface temperature at the previous time step [K]
%t_sfc_p;
% (equal to either T_ice, T_snow or to T_wML)

%  Output (procedure arguments)

% Updated surface temperature [K]
% (equal to the updated value of either T_ice, T_snow or T_wML)

%  Local variables of type LOGICAL
% Switch, true = ice does not exist but should be created
l_ice_create=false;
% Switch, true = there is snow above the ice
l_snow_exists=false;
% Switch, true = snow/ice melting from above takes place
l_ice_meltabove=[];

%  Local variables of type INTEGER
% Loop index
i=0;
%  Local variables of type REAL
% Time derivative of T_mnw [K s^{-1}]
d_t_mnw_dt=0;
% Time derivative of T_ice [K s^{-1}]
d_t_ice_dt=0;             
% Time derivative of T_bot [K s^{-1}]
d_t_bot_d=0;         
% Time derivative of T_B1 [K s^{-1}]
d_t_b1_dt=0;              
% Time derivative of h_snow [m s^{-1}]
d_h_snow_dt=0;            
% Time derivative of h_ice [m s^{-1}]
d_h_ice_dt=0;             
% Time derivative of h_ML [m s^{-1}]
d_h_ml_dt=0;              
% Time derivative of H_B1 [m s^{-1}]
d_h_b1_dt=0;              
% Time derivative of C_T [s^{-1}]
d_c_t_dt=0;

%  Local variables of type REAL
% The mean buoyancy frequency in the thermocline [s^{-1}]
n_t_mean=0;

% The ZM96 equilibrium SBL depth scale [m]
zm_h_scale=0;             
% The equilibrium CBL depth scale [m]
conv_equil_h_scale=0;

%  Local variables of type REAL
% If h_ice<h_ice_threshold, use quasi-equilibrium ice model
h_ice_threshold=0;


% Help storage variable
flk_str_1=0;           
% Help storage variable
flk_str_2=0;           
% Dimensionless ratio, used to store intermediate results
r_h_icesnow=0;         
% Dimensionless ratio, used to store intermediate results
r_rho_c_icesnow=0;     
% Dimensionless ratio, used to store intermediate results
r_ti_icesnow=0;        
% Dimensionless ratio, used to store intermediate results
r_tstar_icesnow=0;


% Switches and parameters that configure FLake bottom heat flux scheme
%lflk_botsed_use   = false; %false = don't use flake heat flux scheme 
%rflk_depth_bs_ref = 10.0;

%==============================================================================
%  Start calculations
%------------------------------------------------------------------------------

%_dm
% Security. Set time-rate-of-change of prognostic variables to zero.
% Set prognostic variables to their values at the previous time step.
% (This is to avoid spurious changes of prognostic variables
% when FLake is used within a 3D model, e.g. to avoid spurious generation of ice
% at the neighbouring lake points as noticed by Burkhardt Rockel.)
%_dm

%depth_bs=10; %[m] just set for testing
%t_bs = 273.15; %[K]


%__________________________________________________________________________
% pipe in external frocing for using external ice cover scheme
%--------------------------------------------------------------------------

extincoef_water_typ = FLAKE.extincoef_water_typ;

depth_w      = FLAKE.depth_w;
t_snow_p_flk = FLAKE.t_snow_n_flk;
t_ice_p_flk  = FLAKE.t_ice_n_flk;
t_wml_p_flk  = FLAKE.t_wml_n_flk;
t_mnw_p_flk  = FLAKE.t_mnw_n_flk;
t_bot_p_flk  = FLAKE.t_bot_n_flk;
t_b1_p_flk   = FLAKE.t_b1_n_flk;
h_snow_p_flk = FLAKE.h_snow_n_flk;
h_ice_p_flk  = FLAKE.h_ice_n_flk;
h_ml_p_flk   = FLAKE.h_ml_n_flk;
c_t_p_flk    = FLAKE.c_t_n_flk;
i_w_flk      = FLAKE.i_w_flk;
i_bot_flk    = FLAKE.i_bot_flk;
l_ice_create = FLAKE.l_ice_create;
lat          = FLAKE.latitude;
par_coriolis = 2.*7.2921e-5.*sin(lat*pi./180);


d_t_mnw_dt   = 0.;
d_t_ice_dt   = 0.;
d_t_bot_dt   = 0.;
d_t_b1_dt    = 0.;
d_h_snow_dt  = 0.;
d_h_ice_dt   = 0.;
d_h_ml_dt    = 0.;
d_h_b1_dt    = 0.;
d_c_t_dt     = 0.;

t_snow_n_flk = t_snow_p_flk;
t_ice_n_flk  = t_ice_p_flk;
t_wml_n_flk  = t_wml_p_flk;
t_mnw_n_flk  = t_mnw_p_flk;
t_bot_n_flk  = t_bot_p_flk;
t_b1_n_flk   = t_b1_p_flk;
h_snow_n_flk = h_snow_p_flk;
h_ice_n_flk  = h_ice_p_flk;
h_ml_n_flk   = h_ml_p_flk;
c_t_n_flk    = c_t_p_flk;

c_i_flk=0;
%------------------------------------------------------------------------------
%  Compute fluxes, using variables from the previous time step.
%------------------------------------------------------------------------------

%_dm
% At this point, the heat and radiation fluxes, namely,
% Q_snow_flk, Q_ice_flk, Q_w_flk,
% I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,
% the mean radiation flux over the mixed layer, I_intm_0_h_flk,
% and the mean radiation flux over the thermocline, I_intm_h_D_flk,
% should be known.
% They are computed within 'flake_interface' (or within the driving model)
% and are available to 'flake_driver'
% through the above variables declared in the MODULE 'flake'.
% In case a lake is ice-covered, Q_w_flk is re-computed below.
%_dm

q_snow_flk     = FLAKE.q_snow_flk;
q_ice_flk      = FLAKE.q_ice_flk;
q_w_flk        = FLAKE.q_w_flk;


i_atm_flk      = FLAKE.i_atm_flk;
i_ice_flk      = FLAKE.i_ice_flk;
i_snow_flk     = FLAKE.i_snow_flk;
i_h_flk        = FLAKE.i_h_flk;
i_intm_0_h_flk = FLAKE.i_intm_0_h_flk;
i_intm_h_d_flk = FLAKE.i_intm_h_d_flk;

u_star_w_flk   = FLAKE.u_star_w_flk;
dmsnowdt_flk   = 0; %!!!don't use FLAKE snow scheme

% Heat flux through the ice-water interface
% Ice exists
% if(h_ice_p_flk>=h_ice_min_flk)
%     % Mixed-layer depth is zero, compute flux
%     if(h_ml_p_flk<=h_ml_min_flk)
%         % Flux with linear T(z)
%         q_w_flk = -tpl_kappa_w.*(t_bot_p_flk-t_wml_p_flk)./depth_w;
%         % d\Phi(0)/d\zeta (thermocline)
%         phi_t_pr0_flk = phi_t_pr0_1.*c_t_p_flk-phi_t_pr0_2;
%         % Account for an increased d\Phi(0)/d\zeta
%         q_w_flk = q_w_flk.*max(phi_t_pr0_flk, 1.);
%         %q_w_flk = FLAKE.q_w_flk;%!!! use external q_w_flk
%     else
%         % Mixed-layer depth is greater than zero, set flux to zero
%         q_w_flk = 0.;
%     end
% end;
% % disp([q_w_flk FLAKE.q_ice_flk])

% A generalized heat flux scale
 q_star_flk = q_w_flk + i_w_flk + i_h_flk - 2..*i_intm_0_h_flk;
 
% %Heat flux through the water-bottom sediment interface
% if(lflk_botsed_use)
%     q_bot_flk = -tpl_kappa_w.*(t_b1_p_flk-t_bot_p_flk)./max(h_b1_p_flk, h_b1_min_flk).*phi_b1_pr0;
% else
%     % The bottom-sediment scheme is not used
%     q_bot_flk = 0.;
% end;
% 
% couple to CryoGrid soil scheme
% set external t_b1_p_flk = temperature of first soil cell
% set external h_b1_p_flk = delta K grid of first soil cell
% tpl_kappa_w = heat conductivity is assumed to be water 
%
%q_bot_flk = - k_b1_p_flk .* (t_b1_p_flk - t_bot_p_flk) ./ h_b1_p_flk; 

%q_bot_flk = 0;
q_bot_flk =  FLAKE.q_bot_flk;


% %------------------------------------------------------------------------------
% %  Check if ice exists or should be created.
% %  If so, compute the thickness and the temperature of ice and snow.
% %------------------------------------------------------------------------------
% % 
% %_dm
% % Notice that a quasi-equilibrium ice-snow model is used
% % to avoid numerical instability when the ice is thin.
% % This is always the case when new ice is created.
% %_dm
% 
% %_dev
% % The dependence of snow density and of snow heat conductivity
% % on the snow thickness is accounted for parametrically.
% % That is, the time derivatives of \rho_S and \kappa_S are neglected.
% % The exception is the equation for the snow thickness
% % in case of snow accumulation and no melting,
% % where d\rho_S/dt is incorporated.
% % Furthermore, some (presumably small) correction terms incorporating
% % the snow density and the snow heat conductivity are dropped out.
% % Those terms may be included as better formulations
% % for \rho_S and \kappa_S are available.
% %_dev
% 
% % Default values
% l_ice_create    = false;
% l_ice_meltabove = false;
% 
% 
% % Ice does not exist
% if(h_ice_p_flk<h_ice_min_flk)
%     
%     l_ice_create = t_wml_p_flk<=(tpl_t_f+c_small_flk) & q_w_flk<0.;
%     % Ice does not exist but should be created
%     if(l_ice_create)
%         d_h_ice_dt = -q_w_flk./tpl_rho_i./tpl_l_f;
%         % Advance h_ice
%         h_ice_n_flk = h_ice_p_flk + d_h_ice_dt.*del_time;
%         % Ice temperature
%         t_ice_n_flk = tpl_t_f + h_ice_n_flk.*q_w_flk./tpl_kappa_i./phi_i_pr0_lin;
%         d_h_snow_dt = dmsnowdt_flk./tpl_rho_s_min;
%         % Advance h_snow
%         if (d_h_snow_dt>0)
%             h_snow_n_flk = h_snow_p_flk + d_h_snow_dt.*del_time;
%             % d\Phi_I(1)/d\zeta_I (ice)
%             phi_i_pr1_flk = phi_i_pr1_lin+ phi_i_ast_mr.*min(1., h_ice_n_flk./h_ice_max);
%             r_h_icesnow = phi_i_pr1_flk./phi_s_pr0_lin.*tpl_kappa_i./flake_snowheatconduct(h_snow_n_flk).* h_snow_n_flk./max(h_ice_n_flk, h_ice_min_flk);
%             % Snow temperature
%             t_snow_n_flk = t_ice_n_flk + r_h_icesnow.*(t_ice_n_flk-tpl_t_f);
%         else
%             h_snow_n_flk=0;
%             r_h_icesnow=0;
%             t_snow_n_flk=t_ice_n_flk;
%         end
%                 
%     end
%     
%     % Ice exists
% else
%     
%     % Check if there is snow above the ice
%     l_snow_exists = h_snow_p_flk>=h_snow_min_flk;
%     
%     % T_sfc = T_f, check for melting from above
%     if(t_snow_p_flk>=(tpl_t_f-c_small_flk))
%         % T_snow = T_ice if snow is absent
%         % There is snow above the ice
%         if(l_snow_exists)
%             % Atmospheric forcing
%             flk_str_1 = q_snow_flk + i_snow_flk - i_ice_flk;
%             % Melting of snow and ice from above
%             if(flk_str_1>=0.)
%                 l_ice_meltabove = true;
%                 d_h_snow_dt =(-flk_str_1./tpl_l_f+dmsnowdt_flk)./flake_snowdensity(h_snow_p_flk);
%                 d_h_ice_dt  = -(i_ice_flk - i_w_flk - q_w_flk)./tpl_l_f./tpl_rho_i;
%             end
%             % No snow above the ice
%         else
%             % Atmospheric forcing + heating from the water
%             flk_str_1 = q_ice_flk + i_ice_flk - i_w_flk - q_w_flk;
%             % Melting of ice from above, snow accumulation may occur
%             if(flk_str_1>=0.)
%                 l_ice_meltabove = true;
%                 d_h_ice_dt  = -flk_str_1./tpl_l_f./tpl_rho_i;
%                 d_h_snow_dt = dmsnowdt_flk./tpl_rho_s_min;
%             end
%         end
%         % Melting from above takes place
%         if(l_ice_meltabove)
%             % Advance h_ice
%             h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt .*del_time;
%             % Advance h_snow
%             h_snow_n_flk = h_snow_p_flk + d_h_snow_dt.*del_time;
%             % Set T_ice to the freezing point
%             t_ice_n_flk  = tpl_t_f;
%             % Set T_snow to the freezing point
%             t_snow_n_flk = tpl_t_f;
%         end
%         
%     end
%     
%     % No melting from above
%     if(~l_ice_meltabove)
%         
%         %d_h_snow_dt = flake_snowdensity(h_snow_p_flk);
%         d_h_snow_dt = 0; %!!! avoid snow cover
%         % Account for d\rho_S/dt
%         if(d_h_snow_dt<tpl_rho_s_max)
%             flk_str_1 = h_snow_p_flk.*tpl_gamma_rho_s./tpl_rho_w_r;
%             flk_str_1 = flk_str_1./(1.-flk_str_1);
%             % Snow density is equal to its maximum value, d\rho_S/dt=0
%         else
%             flk_str_1 = 0.;
%         end
%         % Snow accumulation
%         d_h_snow_dt = dmsnowdt_flk./d_h_snow_dt./(1.+flk_str_1);
%         % Advance h_snow
%         h_snow_n_flk = h_snow_p_flk + d_h_snow_dt.*del_time;
%         
%         % h_ice relative to its maximum value
%         phi_i_pr0_flk = h_ice_p_flk./h_ice_max;
%         % Shape factor (ice)
%         c_i_flk = c_i_lin - c_i_mr.*(1.+phi_i_ast_mr).*phi_i_pr0_flk;
%         % d\Phi_I(1)/d\zeta_I (ice)
%         phi_i_pr1_flk = phi_i_pr1_lin + phi_i_ast_mr.*phi_i_pr0_flk;
%         % d\Phi_I(0)/d\zeta_I (ice)
%         phi_i_pr0_flk = phi_i_pr0_lin - phi_i_pr0_flk;
%         
%         h_ice_threshold = max(1., 2..*c_i_flk.*tpl_c_i.*(tpl_t_f-t_ice_p_flk)./tpl_l_f);
%         h_ice_threshold = phi_i_pr0_flk./c_i_flk.*tpl_kappa_i./tpl_rho_i./tpl_c_i.*h_ice_threshold;
%         % Threshold value of h_ice
%         h_ice_threshold = sqrt(h_ice_threshold.*del_time);
%         h_ice_threshold = min(0.9.*h_ice_max, max(h_ice_threshold, h_ice_min_flk));
%         % h_ice(threshold) < 0.9*H_Ice_max
%         
%         % Use a quasi-equilibrium ice model
%         if(h_ice_p_flk<h_ice_threshold)
%             
%             % Use fluxes at the air-snow interface
%             if(l_snow_exists)
%                 flk_str_1 = q_snow_flk + i_snow_flk - i_w_flk;
%                 % Use fluxes at the air-ice interface
%             else
%                 flk_str_1 = q_ice_flk + i_ice_flk - i_w_flk;
%             end
%             d_h_ice_dt = -(flk_str_1-q_w_flk)./tpl_l_f./tpl_rho_i;
%             % Advance h_ice
%             h_ice_n_flk = h_ice_p_flk + d_h_ice_dt .*del_time;
%             % Ice temperature
%             t_ice_n_flk = tpl_t_f + h_ice_n_flk.*flk_str_1./tpl_kappa_i./phi_i_pr0_flk;
%             % Use a complete ice model
%         else
%             
%             d_h_ice_dt = tpl_kappa_i.*(tpl_t_f-t_ice_p_flk)./h_ice_p_flk.*phi_i_pr0_flk;
%             d_h_ice_dt =(q_w_flk+d_h_ice_dt)./tpl_l_f./tpl_rho_i;
%             
%             % Advance h_ice
%             h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt.*del_time;
%             
%             % Dimensionless parameter
%             r_ti_icesnow = tpl_c_i.*(tpl_t_f-t_ice_p_flk)./tpl_l_f;
%             % Dimensionless parameter
%             r_tstar_icesnow = 1. - c_i_flk;
%             % There is snow above the ice
%             if(l_snow_exists)
%                 r_h_icesnow = phi_i_pr1_flk./phi_s_pr0_lin.*tpl_kappa_i./flake_snowheatconduct(h_snow_p_flk).* h_snow_p_flk./h_ice_p_flk;
%                 r_rho_c_icesnow = flake_snowdensity(h_snow_p_flk).*tpl_c_s./tpl_rho_i./tpl_c_i;
%                 %_dev
%                 %_dm
%                 % These terms should be included as an improved understanding of the snow scheme is gained,
%                 % of the effect of snow density in particular.
%                 %_dm
%                 %_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                           &
%                 %_nu                        + (1.+C_S_lin*h_snow_p_flk/h_ice_p_flk)*R_H_icesnow*R_rho_c_icesnow
%                 %_dev
%                 
%                 % Dimensionless parameter
%                 r_tstar_icesnow = r_tstar_icesnow.*r_ti_icesnow;
%                 
%                 %_dev
%                 %_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                         &
%                 %_nu                        + (1.-R_rho_c_icesnow)*tpl_c_I*T_ice_p_flk/tpl_L_f
%                 %_dev
%                 % Atmospheric fluxes
%                 flk_str_2 = q_snow_flk+i_snow_flk-i_w_flk;
%                 flk_str_1  = c_i_flk.*h_ice_p_flk +(1.+c_s_lin.*r_h_icesnow).*r_rho_c_icesnow.*h_snow_p_flk;
%                 % Effect of snow accumulation
%                 d_t_ice_dt = -(1.-2..*c_s_lin).*r_h_icesnow.*(tpl_t_f-t_ice_p_flk).* tpl_c_s.*dmsnowdt_flk;
%                 % No snow above the ice
%             else
%                 % Dimensionless parameter
%                 r_tstar_icesnow = r_tstar_icesnow.*r_ti_icesnow;
%                 % Atmospheric fluxes
%                 flk_str_2 = q_ice_flk+i_ice_flk-i_w_flk;
%                 flk_str_1  = c_i_flk.*h_ice_p_flk;
%                 d_t_ice_dt = 0.;
%             end
%             % Add flux due to heat conduction
%             d_t_ice_dt = d_t_ice_dt + tpl_kappa_i.*(tpl_t_f-t_ice_p_flk)./h_ice_p_flk.*phi_i_pr0_flk.*(1.-r_tstar_icesnow);
%             % Add flux from water to ice
%             d_t_ice_dt = d_t_ice_dt - r_tstar_icesnow.*q_w_flk;
%             % Add atmospheric fluxes
%             d_t_ice_dt = d_t_ice_dt + flk_str_2;
%             % Total forcing
%             d_t_ice_dt = d_t_ice_dt./tpl_rho_i./tpl_c_i;
%             % dT_ice/dt
%             d_t_ice_dt = d_t_ice_dt./flk_str_1;
%             % Advance T_ice
%             t_ice_n_flk = t_ice_p_flk + d_t_ice_dt.*del_time;
%         end
%         
%         % h_ice relative to its maximum value
%         phi_i_pr1_flk = min(1., h_ice_n_flk./h_ice_max);
%         % d\Phi_I(1)/d\zeta_I (ice)
%         phi_i_pr1_flk = phi_i_pr1_lin + phi_i_ast_mr.*phi_i_pr1_flk;
%         %r_h_icesnow = phi_i_pr1_flk./phi_s_pr0_lin.*tpl_kappa_i./flake_snowheatconduct(h_snow_n_flk).*h_snow_n_flk./max(h_ice_n_flk, h_ice_min_flk);
%         % Snow temperature
%         t_snow_n_flk = t_ice_n_flk;% + r_h_icesnow.*(t_ice_n_flk-tpl_t_f);
%         
%     end
%     
% end
% %t_ice_n_flk=FLAKE.t_ice_cryoGRID;%!!! added to avoid ice scheme
% % Security, limit h_ice by its maximum value
% h_ice_n_flk = min(h_ice_n_flk, h_ice_max);
% 
% % Security, limit the ice and snow temperatures by the freezing point
% t_snow_n_flk = min(t_snow_n_flk, tpl_t_f);
% t_ice_n_flk =  min(t_ice_n_flk,  tpl_t_f);
% 
% %_tmp
% % Security, avoid too low values (these constraints are used for debugging purposes)
% t_snow_n_flk = max(t_snow_n_flk, 73.15);
% t_ice_n_flk =  max(t_ice_n_flk,  73.15);
% %_tmp
% 
% 
% 
% % Remove too thin ice and/or snow
% % Check ice
% if(h_ice_n_flk<h_ice_min_flk)
%     % Ice is too thin, remove it, and
%     h_ice_n_flk = 0.;
%     % set T_ice to the freezing point.
%     t_ice_n_flk = tpl_t_f;
%     % Remove snow when there is no ice, and
%     h_snow_n_flk = 0.;
%     % set T_snow to the freezing point.
%     t_snow_n_flk = tpl_t_f;
%     % 'Exotic' case, ice has been created but proved to be too thin
%     l_ice_create = false;
%     % Ice exists, check snow
% elseif(h_snow_n_flk<h_snow_min_flk) ;
%     % Snow is too thin, remove it,
%     h_snow_n_flk = 0.;
%     % and set the snow temperature equal to the ice temperature.
%     t_snow_n_flk = t_ice_n_flk;
% end

%------------------------------------------------------------------------------
%  Compute the mean temperature of the water column.
%------------------------------------------------------------------------------


%--------------------------------------------------------------------------
%__________________________________________________________________________
% Ice has just been created, set Q_w to zero
% if(l_ice_create)
%     q_w_flk = 0.;
% end;

d_t_mnw_dt =(q_w_flk - q_bot_flk + i_w_flk - i_bot_flk)./tpl_rho_w_r./tpl_c_w./depth_w;
% Advance T_mnw
t_mnw_n_flk = t_mnw_p_flk + d_t_mnw_dt.*del_time;
% Limit T_mnw by the freezing point
t_mnw_n_flk = max(t_mnw_n_flk, tpl_t_f);


%------------------------------------------------------------------------------
%  Compute the mixed-layer depth, the mixed-layer temperature,
%  the bottom temperature and the shape factor
%  with respect to the temperature profile in the thermocline.
%  Different formulations are used, depending on the regime of mixing.
%------------------------------------------------------------------------------

% Ice exists
if(h_ice_n_flk>=h_ice_min_flk)   
    % Limit the mean temperature under the ice by T_r
    t_mnw_n_flk = min(t_mnw_n_flk, tpl_t_r);
    % The mixed-layer temperature is equal to the freezing point
    t_wml_n_flk = tpl_t_f;
     
    % Ice has just been created
    if(l_ice_create)
        % h_ML=D when ice is created      
%        if (h_ml_p_flk>=depth_w-h_ml_min_flk) %!!! changed since a persistent isothermal layer is not observed
            % Set h_ML to zero
            h_ml_n_flk = 0.;
            % Set C_T to its minimum value
            c_t_n_flk = c_t_min;
            % h_ML<D when ice is created
%         else
%             % h_ML remains unchanged
%             h_ml_n_flk = h_ml_p_flk;
%             % C_T (thermocline) remains unchanged
%             c_t_n_flk = c_t_p_flk;
%         end       
        t_bot_n_flk = t_wml_n_flk -(t_wml_n_flk-t_mnw_n_flk)./c_t_n_flk./(1.-h_ml_n_flk./depth_w);
        % Update the bottom temperature                
        
        % Ice exists and T_bot < T_r, molecular heat transfer
    elseif(t_bot_p_flk<tpl_t_r)
                        
        % h_ML remains unchanged
        h_ml_n_flk = h_ml_p_flk;
        % C_T (thermocline) remains unchanged
        c_t_n_flk = c_t_p_flk;
        t_bot_n_flk = t_wml_n_flk -(t_wml_n_flk-t_mnw_n_flk)./c_t_n_flk./(1.-h_ml_n_flk./depth_w);
        %Update the bottom temperature
                                
        % Ice exists and T_bot = T_r, convection due to bottom heating
    else
        % T_bot is equal to the temperature of maximum density
        t_bot_n_flk = tpl_t_r;
               
        if (i_w_flk-q_bot_flk>=0) %!!! i_w_flk-q_bot_flk>=0 (warming under ice)
            % h_ML > 0
            if(h_ml_p_flk>=c_small_flk)
                % C_T (thermocline) remains unchanged
                c_t_n_flk = c_t_p_flk;
                h_ml_n_flk = depth_w.*(1.-(t_wml_n_flk-t_mnw_n_flk)./(t_wml_n_flk-t_bot_n_flk)./c_t_n_flk);
                % Update the mixed-layer depth
                h_ml_n_flk = max(h_ml_n_flk, 0.);
            % h_ML = 0
            else
                % h_ML remains unchanged
                h_ml_n_flk = h_ml_p_flk;
                c_t_n_flk =(t_wml_n_flk-t_mnw_n_flk)./(t_wml_n_flk-t_bot_n_flk);
                % Update the shape factor (thermocline)
                c_t_n_flk = min(c_t_max, max(c_t_n_flk, c_t_min));
            end
            
        else %!!! i_w_flk-q_bot_flk<0 (cooling from the bottom)
            
            %h_ML remains unchanged
            h_ml_n_flk = h_ml_p_flk;
            % C_T (thermocline) remains unchanged
            c_t_n_flk = c_t_p_flk;
            t_bot_n_flk = t_wml_n_flk -(t_wml_n_flk-t_mnw_n_flk)./c_t_n_flk./(1.-h_ml_n_flk./depth_w);
            % Update the bottom temperature
            
        end
                     
    end
    
    % Security, limit the bottom temperature by T_r
    t_bot_n_flk = min(t_bot_n_flk, tpl_t_r);      

% Open water
else
    
    % Generalised buoyancy flux scale and convective velocity scale
    T_water=t_wml_p_flk;
    flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_t_r);
    flk_str_1 = flake_buoypar*q_star_flk/tpl_rho_w_r/tpl_c_w;
    if (flk_str_1<0.)
        % Convection
        w_star_sfc_flk =(-flk_str_1.*h_ml_p_flk).^(1../3.);
    else
        % Neutral or stable stratification
        w_star_sfc_flk = 0.;
    end
    
    %_dm
    % The equilibrium depth of the CBL due to surface cooling with the volumetric heating
    % is not computed as a solution to the transcendental equation.
    % Instead, an algebraic formula is used
    % that interpolates between the two asymptotic limits.
    %_dm
    conv_equil_h_scale = -q_w_flk/max(i_w_flk, c_small_flk);
    % The equilibrium CBL depth scale is only used above T_r
    if(conv_equil_h_scale>0. && conv_equil_h_scale<1.&& t_wml_p_flk>tpl_t_r)
        conv_equil_h_scale = sqrt(6..*conv_equil_h_scale)+ 2..*conv_equil_h_scale./(1.-conv_equil_h_scale);
        conv_equil_h_scale = min(depth_w, conv_equil_h_scale./extincoef_water_typ);
    else
        % Set the equilibrium CBL depth to zero
        conv_equil_h_scale = 0.;
    end;
    
    % Mean buoyancy frequency in the thermocline
    T_water=0.5.*(t_wml_p_flk+t_bot_p_flk);
    flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_t_r);
    n_t_mean = flake_buoypar.*(t_wml_p_flk-t_bot_p_flk);
    if(h_ml_p_flk<=depth_w-h_ml_min_flk) && (n_t_mean>=0) %!!! && (n_t_mean>=0) added to avoid complex numbers
        % Compute N
        n_t_mean = sqrt(n_t_mean./(depth_w-h_ml_p_flk));
    else
        % h_ML=D, set N to zero
        n_t_mean = 0.;
    end
    
    % The rate of change of C_T
    d_c_t_dt = max([w_star_sfc_flk, u_star_min_flk]).^2;
    % Relaxation time scale for C_T
    d_c_t_dt = n_t_mean.*(depth_w-h_ml_p_flk).^2./ c_relax_c./d_c_t_dt;
    % Rate-of-change of C_T
    d_c_t_dt =(c_t_max-c_t_min)./max(d_c_t_dt, c_small_flk);
    
    % Compute the shape factor and the mixed-layer depth,
    % using different formulations for convection and wind mixing
    
    % C_TT, using C_T at the previous time step
    c_tt_flk = c_tt_1.*c_t_p_flk-c_tt_2;
    % C_Q using C_T at the previous time step
    c_q_flk = 2..*c_tt_flk./c_t_p_flk;
    
    
    % Convective mixing
    if(flk_str_1<0.)
        
        % Update C_T, assuming dh_ML/dt>0
        c_t_n_flk = c_t_p_flk + d_c_t_dt.*del_time;
        % Limit C_T
        c_t_n_flk = min(c_t_max, max(c_t_n_flk, c_t_min));
        % Re-compute dC_T/dt
        d_c_t_dt =(c_t_n_flk-c_t_p_flk)./del_time;
        
        % Compute dh_ML/dt
        if(h_ml_p_flk<=depth_w-h_ml_min_flk)
            % Use a reduced entrainment equation (spin-up)
            if(h_ml_p_flk<=h_ml_min_flk)
                d_h_ml_dt = c_cbl_1./c_cbl_2.*max(w_star_sfc_flk, c_small_flk);
                
                %_dbg
                % WRITE*, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', d_h_ML_dt*del_time
                % WRITE*, '         w_*       = ', w_star_sfc_flk
                % WRITE*, '         \beta*Q_* = ', flk_str_1
                %_dbg
                
                % Use a complete entrainment equation
            else
                r_h_icesnow     = depth_w./h_ml_p_flk;
                r_rho_c_icesnow = r_h_icesnow-1.;
                r_ti_icesnow    = c_t_p_flk./c_tt_flk;
                r_tstar_icesnow =(r_ti_icesnow./2.-1.).*r_rho_c_icesnow + 1.;
                d_h_ml_dt = -q_star_flk.*(r_tstar_icesnow.*(1.+c_cbl_1)-1.) - q_bot_flk;
                % Q_* and Q_b flux terms
                d_h_ml_dt = d_h_ml_dt./tpl_rho_w_r./tpl_c_w;
                flk_str_2 =(depth_w-h_ml_p_flk).*(t_wml_p_flk-t_bot_p_flk).*c_tt_2./c_tt_flk.*d_c_t_dt;
                % Add dC_T/dt term
                d_h_ml_dt = d_h_ml_dt + flk_str_2;
                flk_str_2 = i_bot_flk +(r_ti_icesnow-1.).*i_h_flk - r_ti_icesnow.*i_intm_h_d_flk;
                flk_str_2 = flk_str_2 +(r_ti_icesnow-2.).*r_rho_c_icesnow.*(i_h_flk-i_intm_0_h_flk);
                flk_str_2 = flk_str_2./tpl_rho_w_r./tpl_c_w;
                % Add radiation terms
                d_h_ml_dt = d_h_ml_dt + flk_str_2;
                flk_str_2 = -c_cbl_2.*r_tstar_icesnow.*q_star_flk./tpl_rho_w_r./tpl_c_w./max(w_star_sfc_flk, c_small_flk);
                flk_str_2 = flk_str_2 + c_t_p_flk.*(t_wml_p_flk-t_bot_p_flk);
                % dh_ML/dt = r.h.s.
                d_h_ml_dt = d_h_ml_dt./flk_str_2;
            end
            %_dm
            % Notice that dh_ML/dt may appear to be negative
            % (e.g. due to buoyancy loss to bottom sediments and/or
            % the effect of volumetric radiation heating),
            % although a negative generalized buoyancy flux scale indicates
            % that the equilibrium CBL depth has not yet been reached
            % and convective deepening of the mixed layer should take place.
            % Physically, this situation reflects an approximate character of the lake model.
            % Using the self-similar temperature profile in the thermocline,
            % there is always communication between the mixed layer, the thermocline
            % and the lake bottom. As a result, the rate of change of the CBL depth
            % is always dependent on the bottom heat flux and the radiation heating of the thermocline.
            % In reality, convective mixed-layer deepening may be completely decoupled
            % from the processes underneath. In order to account for this fact,
            % the rate of CBL deepening is set to a small value
            % if dh_ML/dt proves to be negative.
            % This is 'double insurance' however,
            % as a negative dh_ML/dt is encountered very rarely.
            %_dm
            
            %_dbg
            % IF(d_h_ML_dt.LT.0.) THEN
            %   WRITE*, 'FLake: negative d_h_ML_dt during convection, = ', d_h_ML_dt
            %   WRITE*, '                d_h_ML_dt*del_time = ', MAX(d_h_ML_dt, c_small_flk)*del_time
            %   WRITE*, '         u_*       = ', u_star_w_flk
            %   WRITE*, '         w_*       = ', w_star_sfc_flk
            %   WRITE*, '         h_CBL_eqi = ', conv_equil_h_scale
            %   WRITE*, '         ZM scale  = ', ZM_h_scale
            %   WRITE*, '        h_ML_p_flk = ', h_ML_p_flk
            % endif
            %   WRITE*, 'FLake: Convection, = ', d_h_ML_dt
            %   WRITE*, '         Q_*       = ', Q_star_flk
            %   WRITE*, '         \beta*Q_* = ', flk_str_1
            %_dbg
            
            d_h_ml_dt = max(d_h_ml_dt, c_small_flk);
            % Update h_ML
            h_ml_n_flk = h_ml_p_flk + d_h_ml_dt.*del_time;
            % Security, limit h_ML
            h_ml_n_flk = max(h_ml_min_flk, min(h_ml_n_flk, depth_w));
            % Mixing down to the lake bottom
        else
            h_ml_n_flk = depth_w;
            %disp('convec 2 bottom')
        end
        
    % Wind mixing
    else
        %disp(' wind mixing')
        % The surface friction velocity
        d_h_ml_dt = max(u_star_w_flk, u_star_min_flk);   
        zm_h_scale =(abs(par_coriolis)./c_sbl_zm_n + n_t_mean./c_sbl_zm_i).*d_h_ml_dt.^2;
        zm_h_scale = zm_h_scale + flk_str_1./c_sbl_zm_s;
        zm_h_scale = max(zm_h_scale, c_small_flk);
        zm_h_scale = d_h_ml_dt.^3./zm_h_scale;
        % The ZM96 SBL depth scale
        zm_h_scale = max(h_ml_min_flk, min(zm_h_scale, h_ml_max_flk));
        % Equilibrium mixed-layer depth
        zm_h_scale = max(zm_h_scale, conv_equil_h_scale);
        
        %_dm
        % In order to avoid numerical discretization problems,
        % an analytical solution to the evolution equation
        % for the wind-mixed layer depth is used.
        % That is, an exponential relaxation formula is applied
        % over the time interval equal to the model time step.
        %_dm
        
        d_h_ml_dt = c_relax_h.*d_h_ml_dt./zm_h_scale.*del_time;
        % Update h_ML
        h_ml_n_flk = zm_h_scale -(zm_h_scale-h_ml_p_flk).*exp(-d_h_ml_dt);
        % Limit h_ML
        h_ml_n_flk = max(h_ml_min_flk, min(h_ml_n_flk, depth_w));
        % Re-compute dh_ML/dt
        d_h_ml_dt =(h_ml_n_flk-h_ml_p_flk)./del_time;
        % Mixed-layer retreat or stationary state, dC_T/dt<0
        if(h_ml_n_flk<=h_ml_p_flk)
            d_c_t_dt = -d_c_t_dt;
        end
        % Update C_T
        c_t_n_flk = c_t_p_flk + d_c_t_dt.*del_time;
        % Limit C_T
        c_t_n_flk = min(c_t_max, max(c_t_n_flk, c_t_min));
        % Re-compute dC_T/dt
        d_c_t_dt =(c_t_n_flk-c_t_p_flk)./del_time;
        
        %_dbg
        % WRITE*, 'FLake: wind mixing: d_h_ML_dt*del_time = ', d_h_ML_dt*del_time
        % WRITE*, '              h_CBL_eqi = ', conv_equil_h_scale
        % WRITE*, '              ZM scale  = ', ZM_h_scale
        % WRITE*, '              w_*       = ', w_star_sfc_flk
        % WRITE*, '              u_*       = ', u_star_w_flk
        % WRITE*, '             h_ML_p_flk = ', h_ML_p_flk
        %_dbg
        
    end
    
    % Compute the time-rate-of-change of the the bottom temperature,
    % depending on the sign of dh_ML/dt
    % Update the bottom temperature and the mixed-layer temperature
    
    % Mixing did not reach the bottom
    if(h_ml_n_flk<=depth_w-h_ml_min_flk)
        
        % Mixed-layer deepening
        if(h_ml_n_flk>h_ml_p_flk)
            r_h_icesnow     = h_ml_p_flk./depth_w;
            r_rho_c_icesnow = 1.-r_h_icesnow;
            r_ti_icesnow    = 0.5.*c_t_p_flk.*r_rho_c_icesnow+c_tt_flk.*(2..*r_h_icesnow-1.);
            r_tstar_icesnow =(0.5+c_tt_flk-c_q_flk)./r_ti_icesnow;
            r_ti_icesnow    =(1.-c_t_p_flk.*r_rho_c_icesnow)./r_ti_icesnow;
            
            d_t_bot_dt =(q_w_flk-q_bot_flk+i_w_flk-i_bot_flk)./tpl_rho_w_r./tpl_c_w;
            d_t_bot_dt = d_t_bot_dt - c_t_p_flk.*(t_wml_p_flk-t_bot_p_flk).*d_h_ml_dt;
            % Q+I fluxes and dh_ML/dt term
            d_t_bot_dt = d_t_bot_dt.*r_tstar_icesnow./depth_w;
            
            flk_str_2 = i_intm_h_d_flk -(1.-c_q_flk).*i_h_flk - c_q_flk.*i_bot_flk;
            flk_str_2 = flk_str_2.*r_ti_icesnow./(depth_w-h_ml_p_flk)./tpl_rho_w_r./tpl_c_w;
            % Add radiation-flux term
            d_t_bot_dt = d_t_bot_dt + flk_str_2;
            
            flk_str_2 =(1.-c_tt_2.*r_ti_icesnow)./c_t_p_flk;
            flk_str_2 = flk_str_2.*(t_wml_p_flk-t_bot_p_flk).*d_c_t_dt;
            % Add dC_T/dt term
            d_t_bot_dt = d_t_bot_dt + flk_str_2;
            
            % Mixed-layer retreat or stationary state
        else
            % dT_bot/dt=0
            d_t_bot_dt = 0.;
        end
        
        % Update T_bot
        t_bot_n_flk = t_bot_p_flk + d_t_bot_dt.*del_time;
        % Security, limit T_bot by the freezing point
        t_bot_n_flk = max(t_bot_n_flk, tpl_t_f);
        T_water=t_mnw_n_flk;
        tpl = tpl_grav*tpl_a_T*(T_water-tpl_t_r);      
        flk_str_2 =(t_bot_n_flk-tpl_t_r)*flake_buoypar;
        % Security, avoid T_r crossover
        if(flk_str_2<0.)
            t_bot_n_flk = tpl_t_r;
        end
        t_wml_n_flk = c_t_n_flk.*(1.-h_ml_n_flk./depth_w);
        t_wml_n_flk =(t_mnw_n_flk-t_bot_n_flk.*t_wml_n_flk)./(1.-t_wml_n_flk);
        % Security, limit T_wML by the freezing point
        t_wml_n_flk = max(t_wml_n_flk, tpl_t_f);
        
        % Mixing down to the lake bottom
    else
        
        h_ml_n_flk = depth_w;
        t_wml_n_flk = t_mnw_n_flk;
        t_bot_n_flk = t_mnw_n_flk;
        c_t_n_flk = c_t_min;
        
    end
    
end


%------------------------------------------------------------------------------
%  Compute the depth of the upper layer of bottom sediments
%  and the temperature at that depth.
%------------------------------------------------------------------------------

% % The bottom-sediment scheme is used
% if(lflk_botsed_use)
%     
%     % No T(z) maximum (no thermal wave)
%     if(h_b1_p_flk>=depth_bs-h_b1_min_flk)
%         % Set H_B1_p to zero
%         h_b1_p_flk = 0.;
%         % Set T_B1_p to the bottom temperature
%         t_b1_p_flk = t_bot_p_flk;
%     end
%     
%     flk_str_1 = 2..*phi_b1_pr0./(1.-c_b1).*tpl_kappa_w./tpl_rho_w_r./tpl_c_w.*del_time;
%     % Threshold value of H_B1
%     h_ice_threshold = sqrt(flk_str_1);
%     % Limit H_B1
%     h_ice_threshold = min(0.9.*depth_bs, h_ice_threshold);
%     flk_str_2 = c_b2./(1.-c_b2).*(t_bs-t_b1_p_flk)./(depth_bs-h_b1_p_flk);
%     
%     % Use a truncated equation for H_B1(t)
%     if(h_b1_p_flk<h_ice_threshold)
%         % Advance H_B1
%         h_b1_n_flk = sqrt(h_b1_p_flk.^2+flk_str_1);
%         % Re-compute dH_B1/dt
%         d_h_b1_dt =(h_b1_n_flk-h_b1_p_flk)./del_time;
%         % Use a full equation for H_B1(t)
%     else
%         flk_str_1 =(q_bot_flk+i_bot_flk)./h_b1_p_flk./tpl_rho_w_r./tpl_c_w;
%         flk_str_1 = flk_str_1 -(1.-c_b1).*(t_bot_n_flk-t_bot_p_flk)./del_time;
%         d_h_b1_dt =(1.-c_b1).*(t_bot_p_flk-t_b1_p_flk)./h_b1_p_flk + c_b1.*flk_str_2;
%         d_h_b1_dt = flk_str_1./d_h_b1_dt;
%         % Advance H_B1
%         h_b1_n_flk = h_b1_p_flk + d_h_b1_dt.*del_time;
%     end
%     d_t_b1_dt = flk_str_2.*d_h_b1_dt;
%     % Advance T_B1
%     t_b1_n_flk = t_b1_p_flk + d_t_b1_dt.*del_time;
%     
%     %_dbg
%     % WRITE*, 'BS module: '
%     % WRITE*, '  Q_bot   = ', Q_bot_flk
%     % WRITE*, '  d_H_B1_dt = ', d_H_B1_dt
%     % WRITE*, '  d_T_B1_dt = ', d_T_B1_dt
%     % WRITE*, '  H_B1    = ', H_B1_n_flk
%     % WRITE*, '    T_bot = ', T_bot_n_flk
%     % WRITE*, '  T_B1    = ', T_B1_n_flk
%     % WRITE*, '    T_bs  = ',  T_bs
%     %_dbg
%     
%     %_nu
%     % Use a very simplistic procedure, where only the upper layer profile is used,
%     % H_B1 is always set to depth_bs, and T_B1 is always set to T_bs.
%     % Then, the time derivatives are zero, and the sign of the bottom heat flux depends on
%     % whether T_bot is smaller or greater than T_bs.
%     % This is, of course, an oversimplified scheme.
%     %_nu  d_H_B1_dt = 0.
%     %_nu  d_T_B1_dt = 0.
%     %_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
%     %_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1
%     %_nu
%     
%     % H_B1 reached depth_bs, or
%     l_snow_exists = h_b1_n_flk>=depth_bs-h_b1_min_flk;
%     % H_B1 decreased to zero, or
%     % h_b1_n_flk<h_b1_min_flk;
%     % there is no T(z) maximum
%     %(t_bot_n_flk-t_b1_n_flk).*(t_bs-t_b1_n_flk)<=0.;
%     if(l_snow_exists)
%         % Set H_B1 to the depth of the thermally active layer
%         h_b1_n_flk = depth_bs;
%         % Set T_B1 to the climatological temperature
%         t_b1_n_flk = t_bs;
%     end
%     
%     % The bottom-sediment scheme is not used
% else
%     
%     % H_B1 is set to a reference value
%     h_b1_n_flk = rflk_depth_bs_ref;
%     % T_B1 is set to the temperature of maximum density
%     t_b1_n_flk = tpl_t_r;
%     
% end


%------------------------------------------------------------------------------
%  Impose additional constraints.
%------------------------------------------------------------------------------

% In case of unstable stratification, force mixing down to the bottom
T_water=t_mnw_n_flk;
flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_t_r);
flk_str_2 =(t_wml_n_flk-t_bot_n_flk)*flake_buoypar;
if(flk_str_2<0.)
    
    %_dbg
    % WRITE*, 'FLake: inverse (unstable) stratification !!! '
    % WRITE*, '       Mixing down to the bottom is forced.'
    % WRITE*, '  T_wML_p, T_wML_n ', T_wML_p_flk-tpl_T_f, T_wML_n_flk-tpl_T_f
    % WRITE*, '  T_mnw_p, T_mnw_n ', T_mnw_p_flk-tpl_T_f, T_mnw_n_flk-tpl_T_f
    % WRITE*, '  T_bot_p, T_bot_n ', T_bot_p_flk-tpl_T_f, T_bot_n_flk-tpl_T_f
    % WRITE*, '  h_ML_p,  h_ML_n  ', h_ML_p_flk,          h_ML_n_flk
    % WRITE*, '  C_T_p,   C_T_n   ', C_T_p_flk,           C_T_n_flk
    %_dbg
    
    h_ml_n_flk = depth_w;
    t_wml_n_flk = t_mnw_n_flk;
    t_bot_n_flk = t_mnw_n_flk;
    c_t_n_flk = c_t_min;
    %disp('unstable')
    
end;


%------------------------------------------------------------------------------
%  Update the surface temperature.
%------------------------------------------------------------------------------

if(h_snow_n_flk>=h_snow_min_flk)
    % Snow exists, use the snow temperature
    t_sfc_n = t_snow_n_flk;
elseif(h_ice_n_flk>=h_ice_min_flk)
    % Ice exists but there is no snow, use the ice temperature
    t_sfc_n = t_ice_n_flk;
else
    % No ice-snow cover, use the mixed-layer temperature
    t_sfc_n = t_wml_n_flk;
end

%------------------------------------------------------------------------------
%  End calculations
%==============================================================================

FLAKE.t_snow_n_flk =    t_snow_n_flk;
FLAKE.t_ice_n_flk  =    t_ice_n_flk;
FLAKE.t_ice_p_flk  =    t_ice_p_flk;
FLAKE.t_wml_n_flk  =    t_wml_n_flk;
FLAKE.t_wml_p_flk  =    t_wml_p_flk;
FLAKE.t_mnw_n_flk  =    t_mnw_n_flk;
FLAKE.t_mnw_p_flk  =    t_mnw_p_flk;
FLAKE.d_t_mnw_dt   =    d_t_mnw_dt;
FLAKE.t_bot_n_flk  =    t_bot_n_flk;
FLAKE.t_bot_p_flk  =    t_bot_p_flk;
FLAKE.h_snow_n_flk =    h_snow_n_flk;
FLAKE.h_ice_n_flk  =    h_ice_n_flk;
FLAKE.h_ml_n_flk   =    h_ml_n_flk;
FLAKE.c_t_n_flk    =    c_t_n_flk;
FLAKE.c_i_flk      =    c_i_flk;
FLAKE.q_bot_flk    =    q_bot_flk;
FLAKE.d_h_ice_dt   =    d_h_ice_dt;
FLAKE.q_w_flk      =    q_w_flk;


