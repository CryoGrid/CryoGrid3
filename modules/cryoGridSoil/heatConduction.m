function SEB = heatConduction(T, k_eff, GRID, PARA, SEB)
              
Q=PARA.soil.Qgeo;                         
cT_delta = GRID.general.cT_delta;
cT_cellAboveSurface = GRID.air.cT_domain_lb;

% Checks 1
assert(sum(imag(k_eff))==0,'ERROR heatConduction : k_eff is complex');
assert(isnan(GRID.air.cT_domain_lb)==0,'ERROR heatConduction : GRID.air.cT_domain_lb; is NaN');
assert(sum(imag(cT_delta))==0,'ERROR heatConduction : cT delta is complex');

dE_dt=T.*0;
dE_dt(cT_cellAboveSurface+1)= k_eff(cT_cellAboveSurface+2).*(T(cT_cellAboveSurface+2)-T(cT_cellAboveSurface+1))./cT_delta(cT_cellAboveSurface+1) ;

% Checks 2
assert(sum(cT_delta(cT_cellAboveSurface+1)==0)==0,'ERROR heatConduction : Nul divisor at flag 1');
assert(sum(imag(dE_dt))==0,'dE_dt is complex at flag 1')

dE_dt(cT_cellAboveSurface+2:end-1)=  (k_eff(cT_cellAboveSurface+3:end-1).*(T(cT_cellAboveSurface+3:end)-T(cT_cellAboveSurface+2:end-1))./cT_delta(cT_cellAboveSurface+2:end) -...
    k_eff(cT_cellAboveSurface+2:end-2).*(T(cT_cellAboveSurface+2:end-1)-T(cT_cellAboveSurface+1:end-2))./cT_delta(cT_cellAboveSurface+1:end-1));

% Checks 3
assert(sum(cT_delta(cT_cellAboveSurface+2:end)==0)==0,'ERROR heatConduction : Nul divisor at flag 2');
assert(sum(cT_delta(cT_cellAboveSurface+1:end-1)==0)==0,'ERROR heatConduction : Nul divisor at flag 3');
assert(sum(imag(dE_dt))==0,'dE_dt is complex at flag 2')

% lower BC (dT_dt=geothermal heat flux)
dE_dt(end)= Q - k_eff(end-1).*(T(end)-T(end-1))./cT_delta(end);

% Checks 4
assert(sum(cT_delta(end)==0)==0,'ERROR heatConduction : Nul divisor at flag 4');
assert(sum(imag(dE_dt))==0,'dE_dt is complex at flag 3')

SEB.dE_dt_cond=dE_dt;