function [z0u, z0t, z0q]=flake_roughnessLength(fetch, u_a, u_star, h_ice_p_flk)

%------------------------------------------------------------------------------
%
% Description:
%
%  Computes the water-surface or the ice-surface roughness lengths
%  with respect to wind velocity, potential temperature and specific humidity.
%
%  The water-surface roughness lengths with respect to wind velocity is computed
%  from the Charnock formula when the surface is aerodynamically rough.
%  A simple empirical formulation is used to account for the dependence
%  of the Charnock parameter on the wind fetch.
%  When the flow is aerodynamically smooth, the roughness length with respect to
%  wind velocity is proportional to the depth of the viscous sub-layer.
%  The water-surface roughness lengths for scalars are computed using the power-law
%  formulations in terms of the roughness Reynolds number (Zilitinkevich et al. 2001).
%  The ice-surface aerodynamic roughness is taken to be constant.
%  The ice-surface roughness lengths for scalars
%  are computed through the power-law formulations
%  in terms of the roughness Reynolds number (Andreas 2002).
%
%
% Current Code Owner: DWD, Dmitrii Mironov
%  Phone:  +49-69-8062 2705
%  Fax:    +49-69-8062 3721
%  E-mail: dmitrii.mironov@dwd.de
%
% History:
% Version    Date       Name
% ---------- ---------- ----
% 1.00       2005/11/17 Dmitrii Mironov
%  Initial release
% 1.01  2014/08/09     Moritz Langer
%  - modifited to work with CryoGrid3
%
% Code Description:
% Language: Fortran 90.
% Software Standards: 'European Standards for Writing and
% Documenting Exchangeable Fortran 90 Code'.
%==============================================================================

% fetch: Typical wind fetch [m]
% u_a:  Wind speed [m s^{-1}]                           
% u_star: Friction velocity in the surface air layer [m s^{-1}]                           
% h_ice: Ice thickness [m]

%Flake constants
%==============================================================================

  c_Karman      = 0.40      ; % The von Karman constant 
  Pr_neutral    = 1.0       ; % Turbulent Prandtl number at neutral static stability
  Sc_neutral    = 1.0       ; % Turbulent Schmidt number at neutral static stability
  

  z0u_ice_rough = 1.0E-03   ; % Aerodynamic roughness of the ice surface [m] (rough flow)
  c_z0u_smooth  = 0.1       ; % Constant in the expression for z0u (smooth flow) 
  c_z0u_rough   = 1.23E-02  ; % The Charnock constant in the expression for z0u (rough flow)
  c_z0u_rough_L = 1.00E-01  ; % An increased Charnock constant (used as the upper limit)
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
  re_z0u_thresh = 0.1;        % Threshold value of the roughness Reynolds number
                              % [value from Zilitinkevich, Grachev, and Fairall (200),
                              % currently not used]
   

%  Thermodynamic parameters
  tpl_grav          = 9.81       ; % Acceleration due to gravity [m s^{-2}] 
  tpsf_R_dryair    = 2.8705E+02  ; % Gas constant for dry air [J kg^{-1} K^{-1}]
  tpsf_R_watvap    = 4.6151E+02  ; % Gas constant for water vapour [J kg^{-1} K^{-1}]
  
  
  tpsf_nu_u_a      = 1.50E-05    ; % Kinematic molecular viscosity of air [m^{2} s^{-1}]

  
%  Security constants
  u_wind_min_sf  = 1.0E-02  ; % Minimum wind speed [m s^{-1}]
  h_Ice_min_flk  = 1.0E-3   ; % Minimum ice thickness [m] 
  
%  Useful constants
%  num_1o3_sf = 1./3.; % 1/3

%Flake parameters
%==============================================================================


%==============================================================================
%  Start calculations
%------------------------------------------------------------------------------

% Water surface
if(h_ice_p_flk < h_Ice_min_flk)
    
    % The Charnock parameter as dependent on dimensionless fetch
    % Inverse dimensionless fetch
    c_z0u_fetch = max(u_a, u_wind_min_sf).^2./tpl_grav./fetch;
    c_z0u_fetch = c_z0u_rough + c_z0u_ftch_f.*c_z0u_fetch.^c_z0u_ftch_ex;
    % Limit Charnock parameter
    c_z0u_fetch = min(c_z0u_fetch, c_z0u_rough_L);
    
    % Threshold value of friction velocity
    %u_star_thresh =(c_z0u_smooth./c_z0u_fetch.*tpl_grav.*tpsf_nu_u_a).^num_1o3_sf;
    
    % Surface Reynolds number and its threshold value
    re_s = u_star.^3./tpsf_nu_u_a./tpl_grav;
    re_s_thresh = c_z0u_smooth./c_z0u_fetch;
    
    % Aerodynamic roughness
    if(re_s <= re_s_thresh)
        % Smooth flow
        z0u = c_z0u_smooth.*tpsf_nu_u_a./u_star;
    else
        % Rough flow
        z0u = c_z0u_fetch.*u_star.*u_star./tpl_grav;
    end
    % Roughness for scalars
    z0q = c_z0u_fetch.*max(re_s, re_s_thresh);
    z0t = c_z0t_rough_1.*z0q.^c_z0t_rough_3 - c_z0t_rough_2;
    z0q = c_z0q_rough_1.*z0q.^c_z0q_rough_3 - c_z0q_rough_2;
    z0t = z0u.*exp(-c_Karman./Pr_neutral.*z0t);
    z0q = z0u.*exp(-c_Karman./Sc_neutral.*z0q);
    
% Ice surface
else
    
    % Threshold value of friction velocity
    %u_star_thresh = c_z0u_smooth.*tpsf_nu_u_a./z0u_ice_rough;
    
    % Aerodynamic roughness
    z0u = max(z0u_ice_rough, c_z0u_smooth.*tpsf_nu_u_a./u_star);
    
    % Roughness Reynolds number
    re_s = max(u_star.*z0u./tpsf_nu_u_a, 1.0E-07 );
    
    % Roughness for scalars
    if(re_s<=re_z0s_ice_t)
        z0t = c_z0t_ice_b0t + c_z0t_ice_b1t.*log(re_s);
        z0t = min(z0t, c_z0t_ice_b0s);
        z0q = c_z0q_ice_b0t + c_z0q_ice_b1t.*log(re_s);
        z0q = min(z0q, c_z0q_ice_b0s);
    else
        z0t = c_z0t_ice_b0r + c_z0t_ice_b1r.*log(re_s) + c_z0t_ice_b2r.*log(re_s).^2;
        z0q = c_z0q_ice_b0r + c_z0q_ice_b1r.*log(re_s) + c_z0q_ice_b2r.*log(re_s).^2;
    end
    z0t = z0u.*exp(z0t);
    z0q = z0u.*exp(z0q);
    
end


%------------------------------------------------------------------------------
%  End calculations
%==============================================================================
