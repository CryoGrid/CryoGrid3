function [T, GRID, PARA, SEB] = CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_temp, timestep)

    if ~isempty(GRID.snow.cT_domain_ub) %snow cover already exitis

        %----------calculate snow surface albedo ------------------------------

        if max(T(GRID.snow.cT_domain))>=0    % melting conditions
            PARA.snow.albedo = PARA.snow.min_albedo + (PARA.snow.albedo - PARA.snow.min_albedo) ...
                            .* exp(-PARA.snow.tau_f .* timestep.*24.*3600 ./ PARA.snow.tau_1);
        else
            PARA.snow.albedo=max(PARA.snow.albedo-PARA.snow.tau_a.*timestep.*24.*3600./PARA.snow.tau_1, PARA.snow.min_albedo);
        end

        %--------SEB.sublimation-----------------------------------------------

        GRID.snow.Snow_i(GRID.snow.cT_domain_ub) = GRID.snow.Snow_i(GRID.snow.cT_domain_ub) ...
                                                 - SEB.Qe.*timestep.*24.*3600./PARA.constants.L_sg./PARA.constants.rho_i;

        nonAirFractionUppermostGridCell = (GRID.snow.Snow_i(GRID.snow.cT_domain_ub)+GRID.snow.Snow_w(GRID.snow.cT_domain_ub))./...
            (GRID.snow.Snow_i(GRID.snow.cT_domain_ub)+GRID.snow.Snow_w(GRID.snow.cT_domain_ub)+GRID.snow.Snow_a(GRID.snow.cT_domain_ub));

        GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = GRID.snow.Snow_a(GRID.snow.cT_domain_ub) ... 
                                                 - ( SEB.Qe.*timestep.*24.*3600./PARA.constants.L_sg./PARA.constants.rho_i./nonAirFractionUppermostGridCell ...
                                                 -  SEB.Qe.*timestep.*24.*3600./PARA.constants.L_sg./PARA.constants.rho_i); 

        %---------- melt and infiltration -------------------------------------
        if max(T(GRID.snow.cT_domain))>0 || FORCING.i.rainfall>0 || sum(GRID.snow.Snow_w)>0  %cases when melt or infiltration occur

            [T(GRID.snow.cT_domain),...
               GRID.snow.Snow_i(GRID.snow.cT_domain),...
               GRID.snow.Snow_w(GRID.snow.cT_domain),...
               GRID.snow.Snow_a(GRID.snow.cT_domain),...
               newMelt] = ...
                snowMelt(T(GRID.snow.cT_domain),...
                         GRID.snow.Snow_i(GRID.snow.cT_domain),...
                         GRID.snow.Snow_w(GRID.snow.cT_domain),...
                         GRID.snow.Snow_a(GRID.snow.cT_domain),...
                         FORCING.i.rainfall.*timestep./1000,...
                         c_temp(GRID.snow.cT_domain),...
                         PARA);

        end

         %-------- add the new snow to the upper most snow cell in the case of a exisiting snow cover -------------------
        if isempty(PARA.snow.maxSnow)
            deltaSnow_i = max(0, FORCING.i.snowfall.*timestep./1000);
        else
            snowHeight = abs( GRID.general.K_grid(GRID.snow.cT_domain_ub) - GRID.general.K_grid(GRID.snow.cT_domain_lb+1) );
            if snowHeight>PARA.snow.maxSnow
                warning(' excess snow occurs ');
            end
            deltaSnow_i = max( [ 0, ...
                               min( [ FORCING.i.snowfall.*timestep./1000, ...
                                      (PARA.snow.maxSnow - snowHeight ) .* PARA.snow.rho_snow ./ PARA.constants.rho_w ] ) ] ); %ensures that no more than maxSnow can accumulate
        end

        GRID.snow.Snow_i(GRID.snow.cT_domain_ub) = GRID.snow.Snow_i(GRID.snow.cT_domain_ub) ...
                                                 + deltaSnow_i;

        GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = GRID.snow.Snow_a(GRID.snow.cT_domain_ub) ... 
                                                + (deltaSnow_i./(PARA.snow.rho_snow./1000) ...
                                                - deltaSnow_i);    


    else %no snow cover
        %---------- add the new snow into initial SWE variable in case of no snow cover------------------
        GRID.snow.SWEinitial = GRID.snow.SWEinitial + FORCING.i.snowfall.*timestep./1000 - GRID.snow.SWEinitial.*0.1.*timestep;
    end

    %--------- update albedo after fresh fallen snow --------------------------
    % determine time of last snowfall for albedo calculation
    SEB.newSnow =  SEB.newSnow-SEB.newSnow.*0.1.*timestep + FORCING.i.snowfall.*timestep./1000;  
    if SEB.newSnow>= PARA.technical.SWEperCell/2
        PARA.snow.albedo=PARA.snow.max_albedo;
        SEB.newSnow=0;
    end
end