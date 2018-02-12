function [FLAKE GRID] = initializeLAKE(GRID, PARA);


%---- flake initialization ------------------------------------------------
FLAKE.t_snow_n_flk=0+273.15;
FLAKE.t_ice_n_flk=0+273.15;
FLAKE.t_wml_n_flk=6+273.15;
FLAKE.t_mnw_n_flk=5+273.15;
FLAKE.t_bot_n_flk=4+273.15;
FLAKE.t_b1_n_flk=7+273.15;
FLAKE.h_snow_n_flk=0;
FLAKE.h_ice_n_flk=0;
FLAKE.h_ml_n_flk=3;
FLAKE.h_b1_n_flk=10;
FLAKE.c_t_n_flk=0;
FLAKE.h_ice_n_flk=0;

FLAKE.i_w_flk=0;
FLAKE.i_bot_flk=0;
FLAKE.q_w_flk=0;
FLAKE.q_bot_flk=0;

FLAKE.fetch=100;
FLAKE.depth_w=PARA.water.depth;

FLAKE.d_h_ice_dt  = 0;
FLAKE.q_ice_water = 0;
FLAKE.extincoef_water_typ=PARA.water.extinction;
FLAKE.latitude=PARA.location.latitude;