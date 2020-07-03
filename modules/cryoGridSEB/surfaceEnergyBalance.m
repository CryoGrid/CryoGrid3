function SEB=surfaceEnergyBalance(T, wc, FORCING, GRID, PARA, SEB);


Lstar=mean(SEB.L_star);

sigma=5.67e-8; %Stefan-Boltzmann const.
z=PARA.technical.z;

Qh=real(Q_h(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, FORCING.i.p, FORCING.i.q));



%______here SW radiation is calculated_____________________________________
dE_dt=GRID.general.cT_grid.*0;
Qsolar=GRID.general.cT_grid.*0;

dE_dt(GRID.air.cT_domain_lb+1,1)=(1-PARA.surf.albedo).*FORCING.i.Sin;
%------ snow surface (solid state green house effect) ---------------------
if ~isempty(GRID.snow.cT_domain_ub)    
    beta=PARA.snow.extinction;        
    Qsolar(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb+1) = dE_dt(GRID.snow.cT_domain_ub) .* exp(-beta.*(GRID.general.K_grid(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb+1)-GRID.general.K_grid(GRID.snow.cT_domain_ub)));
    dE_dt(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb) = -Qsolar(GRID.snow.cT_domain_ub+1:GRID.snow.cT_domain_lb+1) + Qsolar(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb);
    %put the rest to cell below snow
    dE_dt(GRID.snow.cT_domain_lb+1) = Qsolar(GRID.snow.cT_domain_lb+1);
end


%__________________________________________________________________________
Sout = PARA.surf.albedo*FORCING.i.Sin;
Lout = PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4+(1-PARA.surf.epsilon).*FORCING.i.Lin;
Qnet = FORCING.i.Sin-Sout + FORCING.i.Lin - Lout ;
Qg   = Qnet-Qh-Qe;


%calculate ET

if ~isempty(GRID.snow.cT_domain_ub) || T(GRID.soil.cT_domain_ub)<=0   %snow cover or uppermost grid cell frozen
    Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p));
else
    Qe_pot=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, 0, FORCING.i.p));
    fraction
    
end


dE_dt(GRID.air.cT_domain_lb+1) = dE_dt(GRID.air.cT_domain_lb+1) ...
                               + PARA.surf.epsilon.*FORCING.i.Lin ...
                               - PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4 ...
                               - Qh - Qe;  %Qe positive: cooling of soil =>epaporation/subl. => loss of SWE


                           
                           
                           
SEB.dE_dt_SEB = dE_dt;
SEB.Qnet = Qnet;
SEB.Qh = Qh;
SEB.Qe = Qe;
SEB.Qg = Qg;
SEB.Sout = Sout;
SEB.Lout = Lout;

     