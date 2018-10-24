function Q_e = Q_eq(uz, z, z0, q, Tz, T_surf, Lstar, rs, p, PARA)
%q=specific humidity [kg/kg]
Tz=Tz+273.15;
T_surf=T_surf+273.15;
%p=1005*100; %air preassure [Pa] 
%rho=1.293;
%rho = (p-(p.*q))./(287.058.*Tz) + (p.*q)./(461.495.*Tz); %air density [kg m^(-3)]
rho = p./(PARA.constants.R_a*Tz); %air density [kg m^(-3)]
%cp = PARA.constants.c_a / PARA.constants.rho_a; % 1005; needs to be volumetric as density is calculated temperature-dependent
kappa = PARA.constants.kappa; %0.4;
L_w=1000.*(2500.8 - 2.36.*(T_surf-273.15));   % [J/kg] latent heat of evaporation of water % JAN: why here L_sl=L_sl(T) assumed temperature-dependent? Ref for formula? apparaently from Wikipedia
L_i=PARA.constants.L_sg;  %[J/kg] %1e3.*2834.1; %latent heat of sublimation
%g = PARA.constants.g; %9.81;
%sigma = PARA.constants.sigma; %5.67e-8;




if T_surf<=273.15
    Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)).*(q-satPresIce(T_surf)./p)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar));
else
    Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)).*(q-satPresWater(T_surf)./p)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar)  ...
        + rs.*uz.*kappa.^2./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)));
end

%Q_e = -rho.*L.*kappa.*uz.*kappa./(log(z./z0)).*(RH.*satPresIce(Tz)-satPresIce(T_surf))./p./(log(z./z0));