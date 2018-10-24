function [SEB, dwc_dt]=surfaceEnergyBalanceInfiltration(T, wc, FORCING, GRID, PARA, SEB)


Lstar=mean(SEB.L_star);

sigma=PARA.constants.sigma; %5.67e-8; %Stefan-Boltzmann const.
L=PARA.constants.L_lg.*PARA.constants.rho_w;
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

%__________________________________________________________________________
Sout = PARA.surf.albedo*FORCING.i.Sin;
Lout = PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4 + (1-PARA.surf.epsilon).*FORCING.i.Lin;
Qnet = FORCING.i.Sin-Sout + FORCING.i.Lin - Lout ;



%calculate ET
if PARA.modules.infiltration
    
    % snow cover or uppermost grid cell frozen --> no ET ; this includes the case of a frozen water body
    if ~isempty(GRID.snow.cT_domain_ub) || T(GRID.soil.cT_domain_ub)<=0
        Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
        % unfrozen water body at surface
    elseif GRID.lake.unfrozenWaterSurface
        Qe=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, PARA.surf.rs, FORCING.i.p, PARA));
        dwc_dt(1)=-Qe./L; %in m water per sec, this can be evaporation or condensation
        
        % unfrozen soil surface
    else
        Qe_pot=real(Q_eq(FORCING.i.wind, z, PARA.surf.z0, FORCING.i.q, FORCING.i.Tair, T(GRID.air.cT_domain_lb+1), Lstar, 0, FORCING.i.p, PARA));  %potential ET
        if Qe_pot>0
                % determine index of soil cell to which E and T occur
                i_ALD = PARA.location.bottomBucketSoilcTIndex; % this corresponds to the ALD or gives a reasonable maximum
                i_Emax = GRID.soil.E_lb;
                i_Tmax = GRID.soil.T_lb;
                i_E = min( i_Emax, i_ALD );
                i_T = min( i_Tmax, i_ALD );
                r = PARA.soil.ratioET;
                % cell-wise "efficiencies"
                fraction_E = getE_fraction( T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_E-1), wc(1:i_E), PARA.soil.fieldCapacity );
                fraction_T = getT_fraction( T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_T-1), wc(1:i_T), PARA.soil.fieldCapacity );
                % grid cell heights for weighting
                K_delta_E = GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_E-1);
                K_delta_T = GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+i_T-1);
                % total efficiencies
                efficiency_E = sum( fraction_E .* K_delta_E ) ./ sum( K_delta_E );
                efficiency_T = sum( fraction_T .* K_delta_T ) ./ sum( K_delta_T );
                % actual Qe
                Qe = min( [ efficiency_E + efficiency_T * r / (1-r), 1 ] ) * Qe_pot;
                % actual Qe_E and Qe_T partition
                Qe_E = Qe * efficiency_E / ( efficiency_E + efficiency_T * r / (1-r) );
                Qe_T = Qe * efficiency_T * r / (1-r) / ( efficiency_E + efficiency_T * r / (1-r) );
                % associated changes of water amounts in [m/s]
                dwc_dt(1:i_E) = dwc_dt(1:i_E) + -Qe_E ./ L .* fraction_E .* K_delta_E ./ sum( fraction_E .* K_delta_E );
                dwc_dt(1:i_T) = dwc_dt(1:i_T) + -Qe_T ./ L .* fraction_T .* K_delta_T ./ sum( fraction_T .* K_delta_T );          
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
SEB.dE_dt_SEB = dE_dt;
SEB.Qnet = Qnet;
SEB.Qh = Qh;
SEB.Qe = Qe;
SEB.Qg = Qg;
SEB.Sout = Sout;
SEB.Lout = Lout;

