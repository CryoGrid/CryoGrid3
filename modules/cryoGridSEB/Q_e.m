function Q_e = Q_eq(uz, z, z0, RH, Tz, T_surf, Lstar, rs)

rho=1.293;
cp=1005;
L=2.8*10^6;
kappa=0.4;
g=9.81;
p=1005*100;
sigma=5.67e-8;
T_surf=T_surf+273.15;
Tz=Tz+273.15;

if Tz<=273.15
    Q_e = -rho.*L.*kappa.*uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)).*((RH.*satPresIce(Tz)-satPresIce(T_surf))./p)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar));
else
    Q_e = -rho.*L.*kappa.*uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)).*((RH.*satPresWater(Tz)-satPresWater(T_surf))./p)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar)  ...
        + rs.*uz.*kappa.^2./(log(z./z0)- psi_M(z./Lstar, z0./Lstar))  );
end

%Q_e = -rho.*L.*kappa.*uz.*kappa./(log(z./z0)).*(RH.*satPresIce(Tz)-satPresIce(T_surf))./p./(log(z./z0));