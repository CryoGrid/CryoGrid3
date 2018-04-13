function SEB = L_star(FORCING, PARA, SEB)


uz=FORCING.i.wind;
z=PARA.technical.z;
z0=PARA.surf.z0;
Tz=FORCING.i.Tair+273.15;
p=FORCING.i.p;
Qh=SEB.Qh;
Qe=SEB.Qe;
Lstar=mean(SEB.L_star);          

    
rho = p./(PARA.constants.R_a.*Tz); %air density [kg m^(-3)]
cp = PARA.constants.c_a / PARA.constants.rho_a; %1005;
L=1000.*(2500.8 - 2.36.*(Tz-273.15));  %latent heat of evaporation of water, ref for formula? --> apparently from Wikipedia
kappa = PARA.constants.kappa; %0.4;
g = PARA.constants.g; %9.81;

% potential error here, pressure set to 1005 hPa, must be changed? 

u_star = real(uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)));
L_star = real(-rho.*cp.*Tz./kappa./g.*u_star.^3./(Qh + 0.61.*cp./L.*Tz.*Qe));
L_star=(abs(L_star)<1e-6).*L_star./abs(L_star).*1e-6 + (abs(L_star)>=1e-6).*L_star;  % lower limit for Lstar (before 1e-7 and 1e-5)
L_star=(abs(L_star)>1e6).*L_star./abs(L_star).*1e6 + (abs(L_star)<=1e6).*L_star;    % introduced upper limit for Lstar


SEB.u_star = u_star; 
SEB.L_star=[SEB.L_star(2:end) L_star];  