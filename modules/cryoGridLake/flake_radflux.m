
function [i_atm_flk, i_w_flk, i_ice_flk, i_snow_flk, i_h_flk, i_bot_flk, i_intm_0_h_flk, i_intm_h_d_flk]=flake_radflux(flake_SWnet,FLAKE,PARA) %#codegen

%FlakeParameters
% Security constants
  h_ice_min_flk  = 1.0E-3         ; % Minimum ice thickness [m]
  h_ml_min_flk   = 1.0E-2         ; % Minimum mixed-layer depth [m]

%--------------------------------------------------------------------------
blueice_extincoef_optic=PARA.ice.extinction;
blueice_frac_optic=1;

water_extincoef_optic=PARA.water.extinction;
water_frac_optic=1;

% transwater_extincoef_optic=0.3;
% transwater_frac_optic=1;
% 
% drysnow.extincoef_optic=25;
% drysnow.frac_optic=1;
% 
% meltingsnow.extincoef_optic=15;
% meltingsnow.frac_optic=1;
%--------------------------------------------------------------------------

i_atm_flk=flake_SWnet; %net sw radiation at the top of the ice cover
h_ice_p_flk=0; %set to zero due to external ice cover scheme  
h_ml_p_flk=FLAKE.h_ml_n_flk;



if (h_ice_p_flk >= h_ice_min_flk) % ice exists (snow is excluded)
    i_snow_flk = i_atm_flk;
    i_ice_flk  = i_atm_flk;
    i_bot_flk = blueice_frac_optic .* exp(-blueice_extincoef_optic .* h_ice_p_flk);
    i_w_flk      = i_ice_flk .* i_bot_flk;
else                              % no ice cover
    i_snow_flk   = i_atm_flk;  
    i_ice_flk    = i_atm_flk;
    i_w_flk      = i_atm_flk; 
end
 
 if(h_ml_p_flk >= h_ml_min_flk) %there is a mixed layer
    i_bot_flk = water_frac_optic.*exp(-water_extincoef_optic .* h_ml_p_flk);
    i_h_flk = i_w_flk*i_bot_flk;
 else                           % mixed layer depth is less then a minimum value
    i_h_flk = i_w_flk;
 end
 
 i_bot_flk = water_frac_optic .* exp(-water_extincoef_optic * (FLAKE.depth_w));
 i_bot_flk = i_w_flk.*i_bot_flk;

if(h_ml_p_flk >= h_ml_min_flk)  %integral-mean radiation flux over the mixed layer
    i_intm_0_h_flk = water_frac_optic ./ water_extincoef_optic.* (1 - exp(-water_extincoef_optic .* h_ml_p_flk));
    i_intm_0_h_flk = i_w_flk .* i_intm_0_h_flk ./ h_ml_p_flk;
else
    i_intm_0_h_flk = i_h_flk;
end

if(h_ml_p_flk <= FLAKE.depth_w - h_ml_min_flk) % integral-mean radiation flux over the thermocline
    i_intm_h_d_flk = water_frac_optic ./ water_extincoef_optic .* ( exp(-water_extincoef_optic .* h_ml_p_flk) - exp(-water_extincoef_optic .* FLAKE.depth_w) );
    i_intm_h_d_flk = i_w_flk .* i_intm_h_d_flk ./ (FLAKE.depth_w-h_ml_p_flk);
else
    i_intm_h_d_flk = i_h_flk;
end