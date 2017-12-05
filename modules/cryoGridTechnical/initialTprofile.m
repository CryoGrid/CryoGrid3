function T0 = initialTprofile( Q, conductivity, cT_grid, K_grid,  K_frozen, K_thawed, arraySizeT, Temp_ini)
T0=interp1(Temp_ini(:,1),Temp_ini(:,2),cT_grid,'linear');
last_interp=find(~isnan(T0),1,'last');

T_steady=steadyState(T0(last_interp), Q, conductivity(last_interp:end,:), cT_grid(last_interp:end), K_grid(last_interp:end),  K_frozen, K_thawed, arraySizeT);

T0(last_interp:end)=T_steady';

