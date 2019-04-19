function varargout = FlakeConstants
 
  c_Karman      = 0.40      ; % The von Karman constant 
  Pr_neutral    = 1.0       ; % Turbulent Prandtl number at neutral static stability
  Sc_neutral    = 1.0       ; % Turbulent Schmidt number at neutral static stability
  c_MO_u_stab   = 5.0       ; % Constant of the MO theory (wind, stable stratification)
  c_MO_t_stab   = 5.0       ; % Constant of the MO theory (temperature, stable stratification)
  c_MO_q_stab   = 5.0       ; % Constant of the MO theory (humidity, stable stratification)
  c_MO_u_conv   = 15.0      ; % Constant of the MO theory (wind, convection)
  c_MO_t_conv   = 15.0      ; % Constant of the MO theory (temperature, convection)
  c_MO_q_conv   = 15.0      ; % Constant of the MO theory (humidity, convection)
  c_MO_u_exp    = 0.25      ; % Constant of the MO theory (wind, exponent)
  c_MO_t_exp    = 0.5       ; % Constant of the MO theory (temperature, exponent)
  c_MO_q_exp    = 0.5       ; % Constant of the MO theory (humidity, exponent)
  z0u_ice_rough = 1.0E-03   ; % Aerodynamic roughness of the ice surface [m] (rough flow)
  c_z0u_smooth  = 0.1       ; % Constant in the expression for z0u (smooth flow) 
  c_z0u_rough   = 1.23E-02  ; % The Charnock constant in the expression for z0u (rough flow)
  c_z0u_rough_l = 1.00E-01  ; % An increased Charnock constant (used as the upper limit)
  c_z0u_ftch_f  = 0.70      ; % Factor in the expression for fetch-dependent Charnock parameter
  c_z0u_ftch_ex = 0.3333333 ; % Exponent in the expression for fetch-dependent Charnock parameter
  c_z0t_rough_1 = 4.0       ; % Constant in the expression for z0t (factor) 
  c_z0t_rough_2 = 3.2       ; % Constant in the expression for z0t (factor)
  c_z0t_rough_3 = 0.5       ; % Constant in the expression for z0t (exponent) 
  c_z0q_rough_1 = 4.0       ; % Constant in the expression for z0q (factor)
  c_z0q_rough_2 = 4.2       ; % Constant in the expression for z0q (factor)
  c_z0q_rough_3 = 0.5       ; % Constant in the expression for z0q (exponent)
  c_z0t_ice_b0s = 1.250     ; % Constant in the expression for z0t over ice
  c_z0t_ice_b0t = 0.149     ; % Constant in the expression for z0t over ice
  c_z0t_ice_b1t = -0.550    ; % Constant in the expression for z0t over ice
  c_z0t_ice_b0r = 0.317     ; % Constant in the expression for z0t over ice
  c_z0t_ice_b1r = -0.565    ; % Constant in the expression for z0t over ice
  c_z0t_ice_b2r = -0.183    ; % Constant in the expression for z0t over ice
  c_z0q_ice_b0s = 1.610     ; % Constant in the expression for z0q over ice
  c_z0q_ice_b0t = 0.351     ; % Constant in the expression for z0q over ice
  c_z0q_ice_b1t = -0.628    ; % Constant in the expression for z0q over ice
  c_z0q_ice_b0r = 0.396     ; % Constant in the expression for z0q over ice
  c_z0q_ice_b1r = -0.512    ; % Constant in the expression for z0q over ice
  c_z0q_ice_b2r = -0.180    ; % Constant in the expression for z0q over ice
  re_z0s_ice_t  = 2.5       ; % Threshold value of the surface Reynolds number 
                              % used to compute z0t and z0q over ice (Andreas 2002)
 re_z0u_thresh = 0.1;         % Threshold value of the roughness Reynolds number
                              % [value from Zilitinkevich, Grachev, and Fairall (200),
                              % currently not used]

%  Dimensionless constants 
  c_free_conv   = 0.14;       % Constant in the expressions for fluxes in free convection

%  Dimensionless constants 
  c_lwrad_emis  = 0.99;       % Surface emissivity with respect to the long-wave radiation

%  Thermodynamic parameters
  tpsf_C_stefboltz = 5.67E-08    ; % The Stefan-Boltzmann constant [W m^{-2} K^{-4}]
  tpsf_R_dryair    = 2.8705E+02  ; % Gas constant for dry air [J kg^{-1} K^{-1}]
  tpsf_R_watvap    = 4.6151E+02  ; % Gas constant for water vapour [J kg^{-1} K^{-1}]
  tpsf_c_a_p       = 1.005E+03   ; % Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  tpsf_L_evap      = 2.501E+06   ; % Specific heat of evaporation [J kg^{-1}]
  tpsf_nu_u_a      = 1.50E-05    ; % Kinematic molecular viscosity of air [m^{2} s^{-1}]
  tpsf_kappa_t_a   = 2.20E-05    ; % Molecular temperature conductivity of air [m^{2} s^{-1}]
  tpsf_kappa_q_a   = 2.40E-05    ; % Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

%  Derived thermodynamic parameters
  tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap        ; % Ratio of gas constants (Rd/Rv)
  tpsf_alpha_q  = (1.-tpsf_Rd_o_Rv)/tpsf_Rd_o_Rv     ; % Diemsnionless ratio 

%  Thermodynamic parameters
  P_a_ref             = 1.0E+05   ; % Reference pressure [N m^{-2} = kg m^{-1} s^{-2}]


%  Security constants
  u_wind_min_sf  = 1.0E-02  ; % Minimum wind speed [m s^{-1}]
  u_star_min_sf  = 1.0E-04  ; % Minimum value of friction velocity [m s^{-1}]
  c_accur_sf     = 1.0E-07  ; % A small number (accuracy)
  c_small_sf     = 1.0E-04  ; % A small number (used to compute fluxes)

%  Useful constants
  num_1o3_sf = 1./3.; % 1/3
  
  
  % OUTPUT of all variables
  fieldNames=who;
  nFields = length(fieldNames);
  
  for iField = 1:nFields
      assignin('caller',fieldNames{iField},eval([fieldNames{iField}]));
  end
  
  
  
  
