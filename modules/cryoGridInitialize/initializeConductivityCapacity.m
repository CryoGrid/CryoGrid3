function [c_temp, k_temp, k_eff, lwc_temp] = initializeConductivityCapacity(T, wc, GRID, PARA)

    c_temp = zeros(size(GRID.general.cT_grid));
    k_temp = zeros(size(GRID.general.cT_grid));
    k_eff = zeros(size(GRID.general.K_grid));
    lwc_temp = zeros(size(GRID.general.cT_grid));

    %------- unused grid cells --------------------------------
    c_temp(GRID.air.cT_domain) = PARA.constants.c_a;    %set some value e.g. air
    k_temp(GRID.air.cT_domain) = PARA.constants.k_a;    %set some value e.g. air
    lwc_temp(GRID.air.cT_domain) = 0;   

    %------- soil domain --------------------------------------
    [c_temp(GRID.soil.cT_domain),...
     k_temp(GRID.soil.cT_domain),...
     lwc_temp(GRID.soil.cT_domain)] = readThermalParameters(T(GRID.soil.cT_domain), GRID, PARA);

    %adjust for the unfrozen part of the domain
    c_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*c_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* capacityUnfrozen(wc,GRID,PARA);
    k_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*k_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* conductivityUnfrozen(wc,GRID,PARA);
    lwc_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*lwc_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* wc;

    %------- snow domain --------------------------------------
    c_temp(GRID.snow.cT_domain) = cap_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
                                           GRID.snow.Snow_w(GRID.snow.cT_domain),...
                                           GRID.snow.Snow_a(GRID.snow.cT_domain),...
                                           PARA);

    k_temp(GRID.snow.cT_domain) = cond_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
                                            GRID.snow.Snow_w(GRID.snow.cT_domain),...
                                            GRID.snow.Snow_a(GRID.snow.cT_domain));

    lwc_temp(GRID.snow.cT_domain) = GRID.snow.Snow_w(GRID.snow.cT_domain)./GRID.general.K_delta(GRID.snow.cT_domain);
end