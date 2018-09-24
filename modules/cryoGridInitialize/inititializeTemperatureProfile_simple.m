function T = inititializeTemperatureProfile_simple(GRID, PARA, FORCING)

T = zeros(size(GRID.general.cT_grid));

[~, TINI.idx]=min(abs(FORCING.data.t_span-PARA.technical.starttime));
T(GRID.air.cT_domain)   = FORCING.data.Tair(TINI.idx);
T(GRID.snow.cT_domain)   = -4;
T(GRID.soil.cT_domain)  = interp1(PARA.Tinitial(:,1), PARA.Tinitial(:,2), GRID.general.cT_grid(GRID.soil.cT_domain), 'linear');