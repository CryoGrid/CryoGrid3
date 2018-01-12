function [SEB, dwc_dt]=surfaceEnergyBalanceInfiltration(T, wc, FORCING, GRID, PARA, SEB)


Lstar=mean(SEB.L_star);

sigma=PARA.constants.sigma; %5.67e-8; %Stefan-Boltzmann const.
L=PARA.constants.L_lg.*PARA.constants.rho_w;  %2.8*10^6.*1000;   %check this, this seems to be for ice? JAN: yes, should be smaller!
z=PARA.technical.z;

Qh=real(Q_h(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, FORCING.i.p, FORCING.i.q, PARA));

dwc_dt=wc.*0;


%______here SW radiation is calculated_____________________________________
dE_dt=GRID.general.cT_grid.*0;
Qsolar=GRID.general.cT_grid.*0;

dE_dt(GRID.air.cT_domain_lb+1)=(1-PARA.surf.albedo).*FORCING.i.Sin;
%------ snow surface (solid state green house effect) ---------------------
if ~isempty(GRID.snow.cT_domain_ub)    
    beta=PARA.snow.extinction;        
    Qsolar(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb+1) = dE_dt(GRID.snow.cT_domain_ub) .* exp(-beta.*(GRID.general.K_grid(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb+1)-GRID.general.K_grid(GRID.snow.cT_domain_ub)));
    dE_dt(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb) = -Qsolar(GRID.snow.cT_domain_ub+1:GRID.snow.cT_domain_lb+1) + Qsolar(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb);
    %put the rest to cell below snow
    dE_dt(GRID.snow.cT_domain_lb+1) = Qsolar(GRID.snow.cT_domain_lb+1);    
end

% JAN : here a modification for water bodies would be needed (extinction),
% but this would probably make no difference for summer due to mixing,
% maybe different in winter/spring


%__________________________________________________________________________
Sout = PARA.surf.albedo*FORCING.i.Sin;
Lout = PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4 + (1-PARA.surf.epsilon).*FORCING.i.Lin;
Qnet = FORCING.i.Sin-Sout + FORCING.i.Lin - Lout ;



%calculate ET
if PARA.modules.infiltration
    
    %snow cover or uppermost grid cell frozen --> no ET ; JAN: this includes the case of a frozen water body
    if ~isempty(GRID.snow.cT_domain_ub) || T(GRID.soil.cT_domain_ub)<=0   
        Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
    % unfrozen water body at surface
    elseif GRID.lake.unfrozenWaterSurface 
        Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
        dwc_dt(1)=-Qe./L; %in m water per sec, this can be evaporation or condensation
        
    % JAN: this is the "default" case of an unfrozen soil surface    
    else 
        Qe_pot=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, 0, FORCING.i.p, PARA));  %potential ET 
        if Qe_pot>0
            fraction_T=getET_fraction(T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+GRID.soil.T_lb-1), wc(1:GRID.soil.T_lb), PARA.soil.fieldCapacity, PARA.soil.wiltingPoint);
            fraction_E=getET_fraction(T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+GRID.soil.E_lb-1), wc(1:GRID.soil.E_lb), PARA.soil.fieldCapacity, PARA.soil.residualWC);
            fraction_ET = fraction_T.*PARA.soil.ratioET;
            fraction_ET(1:GRID.soil.E_lb) = fraction_ET(1:GRID.soil.E_lb) + fraction_E.*(1-PARA.soil.ratioET);

            Qe=sum(fraction_ET.*GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+GRID.soil.T_lb-1))./sum(GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+GRID.soil.T_lb-1)).*Qe_pot;
            fraction_ET=fraction_ET.*GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+GRID.soil.T_lb-1)./sum(fraction_ET.*GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+GRID.soil.T_lb-1));
            % sum(fraction_ET) is always 1 
            dwc_dt(1:GRID.soil.T_lb)=-Qe./L.*fraction_ET;    %in m water per sec
        else  %condensation
            Qe=Qe_pot;
            dwc_dt(1)=-Qe./L; %in m water per sec, put everything in uppermost grid cell
        end
    end
else % this is identical to case with snow cover or frozen ground
    Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
end
%ground heat flux
Qg   = Qnet-Qh-Qe;

%surface heat flux (into upper cell, ground heat flux regards also other
%grid cells, should be identical if no snow cover and no evapotranspiration
%occur
dE_dt(GRID.air.cT_domain_lb+1) = dE_dt(GRID.air.cT_domain_lb+1) ...
                               + PARA.surf.epsilon.*FORCING.i.Lin ...
                               - PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4 ...
                               - Qh - Qe;  % Qe positive: cooling of soil => evaporation/subl. => loss of SWE

                           
                          
% fluxes are in [ W / m^2 ]                       
SEB.Qsurf = dE_dt(GRID.air.cT_domain_lb+1);   

% if abs( SEB.Qsurf-Qg ) > 1e-6
%     warning ( ' Qsurf != Qg ' );
% end
                           
SEB.dE_dt_SEB = dE_dt;
SEB.Qnet = Qnet;
SEB.Qh = Qh;
SEB.Qe = Qe;
SEB.Qg = Qg;
SEB.Sout = Sout;
SEB.Lout = Lout;

     