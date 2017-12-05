    TEMPORARY.timestep_sum=TEMPORARY.timestep_sum+(timestep*24*3600)*timestep;
    TEMPORARY.T_sum=TEMPORARY.T_sum+T.*timestep;
    TEMPORARY.Qe_sum=TEMPORARY.Qe_sum+SEB.Qe.*timestep;
	TEMPORARY.Qh_sum=TEMPORARY.Qh_sum+SEB.Qh.*timestep;
	TEMPORARY.Qnet_sum=TEMPORARY.Qnet_sum+SEB.Qnet.*timestep;
    TEMPORARY.Qg_sum=TEMPORARY.Qg_sum+SEB.Qg.*timestep;
    
    %----store in output table --------------------------------------------
    if  t==TEMPORARY.outputTime
	       
        TEMPORARY.counter = TEMPORARY.counter+1;
       
        %average over timespan:
        TEMPORARY.dt_out=t-TEMPORARY.t_last;
        TEMPORARY.t_last=t; 
        
        TEMPORARY.T_out=TEMPORARY.T_sum./TEMPORARY.dt_out;
        
        TEMPORARY.Qh=TEMPORARY.Qh_sum./TEMPORARY.dt_out;
        TEMPORARY.Qe=TEMPORARY.Qe_sum./TEMPORARY.dt_out;
        TEMPORARY.Qnet=TEMPORARY.Qnet_sum./TEMPORARY.dt_out;
        TEMPORARY.Qg=TEMPORARY.Qg_sum./TEMPORARY.dt_out;
        
        TEMPORARY.timestep_out=TEMPORARY.timestep_sum./TEMPORARY.dt_out;             
        
        %reset sum variables
        TEMPORARY.T_sum(:)=0;
        
        TEMPORARY.Qh_sum=0;
        TEMPORARY.Qe_sum=0;
        TEMPORARY.Qnet_sum=0;
        TEMPORARY.Qg_sum=0;     
        
        TEMPORARY.timestep_sum=0;
        
        %store new values in OUT struct -----------------------------------
        OUT.cryoGrid3=[OUT.cryoGrid3 [NaN(GRID.air.cT_domain_lb,1); TEMPORARY.T_out(GRID.air.cT_domain_lb+1:end,1)]];
        OUT.water=[OUT.water [ NaN( GRID.soil.cT_domain_ub-1,1) ; wc ] ];
        OUT.liquidWater=[OUT.liquidWater [ NaN( GRID.soil.cT_domain_ub-1,1) ; lwc_cTgrid(GRID.soil.cT_domain) ] ];
        OUT.TIMESTEP=[OUT.TIMESTEP; TEMPORARY.timestep_out];
        OUT.timestamp=[OUT.timestamp; t]; 
        
        OUT.snow.outSnow_i=[OUT.snow.outSnow_i GRID.snow.Snow_i];
        OUT.snow.outSnow_a=[OUT.snow.outSnow_a GRID.snow.Snow_a];
        OUT.snow.outSnow_w=[OUT.snow.outSnow_w GRID.snow.Snow_w];
               
        % surface energy balance      
        OUT.SEB.Lsta=[OUT.SEB.Lsta; mean(SEB.L_star)];
        OUT.SEB.QE=[OUT.SEB.QE; TEMPORARY.Qe];
        OUT.SEB.QH=[OUT.SEB.QH; TEMPORARY.Qh];
        OUT.SEB.QG=[OUT.SEB.QG; TEMPORARY.Qg];
        OUT.SEB.QNET=[OUT.SEB.QNET; TEMPORARY.Qnet];
        OUT.SEB.Tsurf=[OUT.SEB.Tsurf; TEMPORARY.T_out(GRID.air.cT_domain_lb+1)];
        OUT.SEB.albedo_stored=[OUT.SEB.albedo_stored; PARA.surf.albedo];
                
        % soil     
        OUT.soil.soil{1, size(OUT.soil.soil,2)+1}=[GRID.soil.cT_water GRID.soil.cT_mineral GRID.soil.cT_organic];
        OUT.soil.topPosition=[OUT.soil.topPosition; -GRID.general.K_grid(GRID.soil.cT_domain_ub)];
        
        % snow      
        if ~isempty(GRID.snow.cT_domain_ub)
            TEMPORARY.topPosition=-GRID.general.K_grid(GRID.snow.cT_domain_ub);
            TEMPORARY.botPosition=-GRID.general.K_grid(GRID.snow.cT_domain_lb+1);
        else
            TEMPORARY.topPosition=NaN;
            TEMPORARY.botPosition=NaN;
        end
        OUT.snow.topPosition=[OUT.snow.topPosition; TEMPORARY.topPosition];
        OUT.snow.botPosition=[OUT.snow.botPosition; TEMPORARY.botPosition];
        
        %------------------------------------------------------------------     
        disp([datestr(now,'yyyy-mm-dd HH:MM:SS'),':  at ',datestr(t), ',  Average timestep: ',  num2str(TEMPORARY.timestep_out), ' seconds'])
      
        TEMPORARY.outputTime=round((TEMPORARY.outputTime+PARA.technical.outputTimestep)./PARA.technical.outputTimestep).*PARA.technical.outputTimestep;
     
        %write output files      
        if  round((t-TEMPORARY.saveTime).*48)==0   
            save([run_number '/' run_number '_output' datestr(t,'yyyy')  '.mat'], 'OUT') 
            OUT = generateOUT();  
            TEMPORARY.saveTime=datenum(str2num(datestr(t,'yyyy'))+1, str2num(datestr(t,'mm')), str2num(datestr(t,'dd')), str2num(datestr(t,'HH')), str2num(datestr(t,'MM')), 0);
        end            
    end