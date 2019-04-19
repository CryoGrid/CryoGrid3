function varargout = FlakeParameters

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
  c_t_min       = 0.5             ; % Minimum value of the shape factor C_T (thermocline)
  c_t_max       = 0.8             ; % Maximum value of the shape factor C_T (thermocline)
  phi_t_pr0_1   = 40./3.   ; % Constant in the expression for the T shape-function derivative 
  phi_t_pr0_2   = 20./3.   ; % Constant in the expression for the T shape-function derivative 
  c_tt_1        = 11./18.  ; % Constant in the expression for C_TT (thermocline)
  c_tt_2        = 7./45.   ; % Constant in the expression for C_TT (thermocline)
  c_b1          = 2./3.    ; % Shape factor (upper layer of bottom sediments)
  c_b2          = 3./5.    ; % Shape factor (lower layer of bottom sediments)
  phi_b1_pr0    = 2.              ; % B1 shape-function derivative 
  c_s_lin       = 0.5             ; % Shape factor (linear temperature profile in the snow layer)
  phi_s_pr0_lin = 1.              ; % S shape-function derivative (linear profile) 
  c_i_lin       = 0.5             ; % Shape factor (linear temperature profile in the ice layer)
  phi_i_pr0_lin = 1.              ; % I shape-function derivative (linear profile) 
  phi_i_pr1_lin = 1.              ; % I shape-function derivative (linear profile) 
  phi_i_ast_mr  = 2.              ; % Constant in the MR2004 expression for I shape factor
  c_i_mr        = 1./12.   ; % Constant in the MR2004 expression for I shape factor
  h_ice_max     = 3.       ;          % Maximum ice tickness in 
                                      % the Mironov and Ritter (2004, MR2004) ice model [m] 

%  Security constants
  h_snow_min_flk = 1.0E-5         ; % Minimum snow thickness [m]
  h_ice_min_flk  = 1.0E-3         ; % Minimum ice thickness [m]
  h_ml_min_flk   = 1.0E-2         ; % Minimum mixed-layer depth [m]
  h_ml_max_flk   = 1.0E+3         ; % Maximum mixed-layer depth [m]
  h_b1_min_flk   = 1.0E-3         ; % Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6         ; % Minimum value of the surface friction velocity [m s^{-1}]

%  Security constant(s)
  c_small_flk    = 1.0E-10;            % A small number

%  Thermodynamic parameters
  tpl_grav          = 9.81       ; % Acceleration due to gravity [m s^{-2}]
  tpl_t_r           = 277.13     ; % Temperature of maximum density of fresh water [K]
  tpl_t_f           = 273.15     ; % Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05 ; % Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03    ; % Maximum density of fresh water [kg m^{-3}]
  tpl_rho_i         = 9.1E+02    ; % Density of ice [kg m^{-3}]
  tpl_rho_s_min     = 1.0E+02    ; % Minimum snow density [kg m^{-3}]
  tpl_rho_s_max     = 4.0E+02    ; % Maximum snow density [kg m^{-3}]
  tpl_gamma_rho_s   = 2.0E+02    ; % Empirical parameter [kg m^{-4}]  
                                            % in the expression for the snow density 
  tpl_l_f           = 3.3E+05    ; % Latent heat of fusion [J kg^{-1}]
  tpl_c_w           = 4.2E+03    ; % Specific heat of water [J kg^{-1} K^{-1}]
  tpl_c_i           = 2.1E+03    ; % Specific heat of ice [J kg^{-1} K^{-1}]
  tpl_c_s           = 2.1E+03    ; % Specific heat of snow [J kg^{-1} K^{-1}]
  tpl_kappa_w       = 5.46E-01   ; % Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_i       = 2.29       ; % Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_s_min   = 0.2        ; % Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_s_max   = 1.5        ; % Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_gamma_kappa_s = 1.3;         % Empirical parameter [J m^{-2} s^{-1} K^{-1}] 
                                   % in the expression for the snow heat conductivity 

%==============================================================================

  % OUTPUT of all variables
  fieldNames=who;
  nFields = length(fieldNames);
  
  for iField = 1:nFields
      assignin('caller',fieldNames{iField},eval([fieldNames{iField}]));
  end


