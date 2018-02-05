function PARA = loadConstants( PARA )
    % important natural constants, given in SI units
    PARA.constants.kappa = 0.4;                                     % von Kármán constant [-]
    PARA.constants.sigma = 5.6704e-8;                               % Stefan-Boltzmann constant [ W / (m^2 K^4) ]
    PARA.constants.g = 9.81;                                        % gravitational acceleration [m/s^2]
    PARA.constants.p_0 = 100500;                                    % normal pressure (sea level) [Pa=kg/(m s^2)]
    %water
    PARA.constants.rho_w = 1000;                                    % density of liquid water (and ice) [kg/m^3]
    PARA.constants.c_w = 4200 * PARA.constants.rho_w;               % volumetric heat capacity of water [ J / (m^3 K) ]
    PARA.constants.k_w = 0.57;                                      % heat conductivity of water [ W/(mK) ] [Hillel(1982)]
    %ice
    PARA.constants.rho_i = 1000;                                    % density of ice, assumed to be equal to that of water [kg/m^3]
    PARA.constants.c_i = 1900 * PARA.constants.rho_i;               % volumetric heat capacity of ice [ J / (m^3 K) ]
    PARA.constants.k_i = 2.2;                                       % heat conductivity of ice [ W/(mK) ] [Hillel(1982)]
    %latent heat of water
    PARA.constants.T_f = 273.15;                                    % freezing point of water / zero degree Celsius [K]
    PARA.constants.L_sl = 334e3;                                    % specific latent heat of fusion of water [J/kg]            [AMS]
    PARA.constants.L_lg = 2501e3;                                   % specific latent heat of vaporization of water [J/kg]      [AMS]
    PARA.constants.L_sg = PARA.constants.L_sl + PARA.constants.L_lg;% specific latent heat of sublimation of water [J/kg]       
    %air
    PARA.constants.rho_a = 1.293;                                   % density of air [kg/m^3] @ 0°C, sea level
    PARA.constants.c_a = 1005.7 * PARA.constants.rho_a;             %c_a= 0.00125*10^6;%[J/m^3 K]   % volumetric heat capacity of dry air [J/(m^3 K)] @ 0°C, sea level, isobar
    PARA.constants.k_a = 0.0243;                                    %ka=0.025; [Hillel(1982)]       % heat conductivity of air [ W/(mK)] @ 0 °C, sea level pressure     
    PARA.constants.R_a = 287.058;                                   % specific gas constant of air [ J/(kg K) ]
    % organic
    %PARA.constants.rho_o = 1; % n.a.
    PARA.constants.c_o = 2.5e6; %[J/(K m^3)]                        % volumetric heat capacity of organic material [J/(K m^3)]
    PARA.constants.k_o = 0.25;                                      % heat conductivity of organic material [ W/(mK) ] [Hillel(1982)]
    % mineral
    %PARA.constants.rho_m = 1; % n.a.
    PARA.constants.c_m = 2.0e6; %[J/(K m^3)]                          % volumetric heat capacity of minearal material [J/(K m^3)]  
    PARA.constants.k_m = PARA.soil.kh_bedrock;                      % heat conductivity of mineral material / bedrock [ W/(mK)] (specified above) %km=3.8 %mineral [Hillel(1982)]    
end
