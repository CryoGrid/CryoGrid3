function [c_temp, k_temp, k_eff] = getThermalProperties(T, c_temp, k_temp, k_eff, GRID, PARA)

%c_old = c_temp;   
    %------- unused grid cells --------------------------------------------
    c_temp(GRID.air.cT_domain) = 4e5;   %set some value e.g. air
    k_temp(GRID.air.cT_domain) = 0.025; %set some value e.g. air
       
    
    %internal energy due to changed heat capacity
   %  T = T .* c_old./c_temp;
   %  T_old = T;          
    %------- soil domain --------------------------------------------------
    [c_temp(GRID.soil.cT_domain),...
     k_temp(GRID.soil.cT_domain)] = readThermalParameters(T(GRID.soil.cT_domain), GRID, PARA);   
   
    %------- snow domain --------------------------------------------------
    c_temp(GRID.snow.cT_domain) = cap_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
                                           GRID.snow.Snow_w(GRID.snow.cT_domain),...
                                           GRID.snow.Snow_a(GRID.snow.cT_domain));
   
    k_temp(GRID.snow.cT_domain) = cond_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
                                            GRID.snow.Snow_w(GRID.snow.cT_domain),...
                                            GRID.snow.Snow_a(GRID.snow.cT_domain));                             
                                         
    %------- interpolate conductivity to K-grid ---------------------------                                                                    
    k_eff(2:end-1) = GRID.general.K_delta(1:end-1)./(2.*GRID.general.cT_delta) .* (1./k_temp(1:end-1)).^2 ...
                   + GRID.general.K_delta(2:end)  ./(2.*GRID.general.cT_delta) .* (1./k_temp(2:end)).^2;
             
    k_eff(2:end-1) = k_eff(2:end-1).^(-0.5);      

    k_eff(1)     = k_temp(1);
    k_eff(end)   = k_temp(end);
    
    %------ correct upper most value below air-domain ---------------------
    k_eff(GRID.air.K_domain_lb+1) = k_temp(GRID.air.K_domain_lb+1);