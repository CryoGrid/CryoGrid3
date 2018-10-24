function [c_temp, k_temp, k_eff, lwc_temp] = getThermalPropertiesInfiltration(T, wc, c_temp, k_temp, k_eff, lwc_temp, GRID, PARA)

    %------- unused grid cells --------------------------------------------
    c_temp(GRID.air.cT_domain) = PARA.constants.c_a;
    k_temp(GRID.air.cT_domain) = PARA.constants.k_a;
    lwc_temp(GRID.air.cT_domain) = 0;   
    
    %------- soil domain --------------------------------------------------
    [c_temp(GRID.soil.cT_domain),...
     k_temp(GRID.soil.cT_domain),...
     lwc_temp(GRID.soil.cT_domain)] = readThermalParameters(T(GRID.soil.cT_domain), GRID, PARA);   
   
    %adjust for the unfrozen part of the domain
    c_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*c_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* capacityUnfrozen(wc,GRID,PARA);
    k_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*k_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* conductivityUnfrozen(wc,GRID,PARA);
    lwc_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*lwc_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* wc;
     
    %-------- set higher conductivity for free water ----------------------
    % now done in conductivityUnfrozen
    
    
    %------- snow domain --------------------------------------------------
    c_temp(GRID.snow.cT_domain) = cap_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
                                           GRID.snow.Snow_w(GRID.snow.cT_domain),...
                                           GRID.snow.Snow_a(GRID.snow.cT_domain),...
                                           PARA);
   
    k_temp(GRID.snow.cT_domain) = cond_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
                                            GRID.snow.Snow_w(GRID.snow.cT_domain),...
                                            GRID.snow.Snow_a(GRID.snow.cT_domain));  
    
    lwc_temp(GRID.snow.cT_domain) = GRID.snow.Snow_w(GRID.snow.cT_domain)./GRID.general.K_delta(GRID.snow.cT_domain);
                                         
    %------- interpolate conductivity to K-grid ---------------------------                                                                    
    k_eff(2:end-1) = GRID.general.K_delta(1:end-1)./(2.*GRID.general.cT_delta) .* (1./k_temp(1:end-1)).^2 ...
                   + GRID.general.K_delta(2:end)  ./(2.*GRID.general.cT_delta) .* (1./k_temp(2:end)).^2;
             
    k_eff(2:end-1) = k_eff(2:end-1).^(-0.5);      

    k_eff(1)     = k_temp(1);
    k_eff(end)   = k_temp(end);
    
    %------ correct upper most value below air-domain ---------------------
    k_eff(GRID.air.K_domain_lb+1) = k_temp(GRID.air.K_domain_lb+1);