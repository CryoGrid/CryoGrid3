function Q_h = Q_h(uz, z, z0, Tz, T_surf, Lstar, p, q, PARA)
%q=specific humidity [kg/kg]
Tz=Tz+273.15;
T_surf=T_surf+273.15;
%p=1005*100; %air preassure [Pa] 
%rho=1.293;
%rho = (p-(p.*q))./(287.058.*Tz) + (p.*q)./(461.495.*Tz); %air density [kg m^(-3)]
rho = p./(PARA.constants.R_a*Tz); %air density [kg m^(-3)]
cp = PARA.constants.c_a / PARA.constants.rho_a; % in [ J/(K kg)] 1005; needs to be volumetric as density is calculated temperature-dependent
kappa = PARA.constants.kappa; %0.4;
%g = PARA.constants.g; %9.81;
%sigma = PARA.constants.sigma; %5.67e-8;

Q_h = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)) .* (Tz-T_surf)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar));