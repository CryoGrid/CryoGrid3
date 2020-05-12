function SEB = heatConduction(T, k_eff, GRID, PARA, SEB)
              
Q=PARA.soil.Qgeo;                         
cT_delta = GRID.general.cT_delta;
cT_cellAboveSurface = GRID.air.cT_domain_lb;

dE_dt=T.*0;

dE_dt(cT_cellAboveSurface+1)= k_eff(cT_cellAboveSurface+2).*(T(cT_cellAboveSurface+2)-T(cT_cellAboveSurface+1))./cT_delta(cT_cellAboveSurface+1) ;

dE_dt(cT_cellAboveSurface+2:end-1)=  (k_eff(cT_cellAboveSurface+3:end-1).*(T(cT_cellAboveSurface+3:end)-T(cT_cellAboveSurface+2:end-1))./cT_delta(cT_cellAboveSurface+2:end) -...
    k_eff(cT_cellAboveSurface+2:end-2).*(T(cT_cellAboveSurface+2:end-1)-T(cT_cellAboveSurface+1:end-2))./cT_delta(cT_cellAboveSurface+1:end-1));


% lower BC (dT_dt=geothermal heat flux)
dE_dt(end)= Q - k_eff(end-1).*(T(end)-T(end-1))./cT_delta(end);

SEB.dE_dt_cond=dE_dt;
