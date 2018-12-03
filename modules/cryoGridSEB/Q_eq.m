function Q_e = Q_eq(uz, z, z0, q, Tz, T_surf, Lstar, rs, p, PARA)

Tz=Tz+273.15;
T_surf=T_surf+273.15;

rho = p./(PARA.constants.R_a*Tz); %air density [kg m^(-3)]
kappa = PARA.constants.kappa; %0.4;
L_w=1000.*(2500.8 - 2.36.*(T_surf-273.15));   % [J/kg] latent heat of evaporation of water
L_i=PARA.constants.L_sg;  %[J/kg] %1e3.*2834.1; %latent heat of sublimation

if isnan(psi_M(z./Lstar, z0./Lstar))==1;
    if isinf(Lstar)==1;
        fprintf('Q_eq : psi_M gives NaN, and Lstar is inf\n')
    elseif Lstar==0;
        fprintf('Q_eq : psi_M gives NaN, and Lstar is null\n')
    else
        fprintf('Q_eq : psi_M gives NaN, but Lstar not inf, not null\n')
    end
elseif isnan(psi_H(z./Lstar, z0./Lstar))==1;
    fprintf('Q_eq : psi_H gives NaN\n')
elseif (log(z./z0) - psi_M(z./Lstar, z0./Lstar))==0
    fprintf('Q_eq : devide by 0 psi_M\n')
elseif (log(z./z0)- psi_H(z./Lstar, z0./Lstar))==0;
    fprintf('Q_eq : devide by 0 psi_H\n')
end

if T_surf<=273.15
    Q_e = -rho.*L_i.*kappa.*uz.*kappa  ./  (log(z./z0) - psi_M(z./Lstar, z0./Lstar))  .*  (q-satPresIce(T_surf)./p)   ./  (log(z./z0)- psi_H(z./Lstar, z0./Lstar));
else
    Q_e = -rho.*L_w.*kappa.*uz.*kappa  ./  (log(z./z0) - psi_M(z./Lstar, z0./Lstar))  .*  (q-satPresWater(T_surf)./p) ./  (log(z./z0)- psi_H(z./Lstar, z0./Lstar)  ...
        + rs.*uz.*kappa.^2./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)));
end