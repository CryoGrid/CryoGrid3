function [T, snow_i, snow_w, snow_a, runoff] = snowMelt(T, snow_i, snow_w, snow_a, water_flux, c_temp, PARA)

%------------melt the snow for cells with T>0------------
    runoff=0;

%energyC=sum(T.*c_temp.*(snow_i+snow_w+snow_a))+(sum(snow_w)+water_flux).*3.34e8;

    poreSpace=(snow_w+snow_a)./(snow_i+snow_w+snow_a);
    poreSpace(isnan(poreSpace(:,1)),1)=300/1000;                %JAN: what does this number mean??

    [T, snow_i, snow_w, snow_a]=melt(T, snow_i, snow_w, snow_a, poreSpace, c_temp);
     
 %   energyC2=sum(T.*c_temp.*(snow_i+snow_w+snow_a))+(sum(snow_w)+water_flux).*3.34e8;
    pS=(snow_w+snow_a)./(snow_i+snow_w+snow_a);
    pS(isnan(pS(:,1)),1)=300/1000;

   




%-----------calculate how much water (in m) a snow cell can hold-----
      
    maxWater=maxLiqWater(T(pS>0), snow_i(pS>0), snow_w(pS>0), snow_a(pS>0), poreSpace(pS>0), c_temp(pS>0));
    
 
        
    if water_flux>0 || (~isempty(maxWater) && min(maxWater)<0)   %infiltration occurs
      
%------------infiltrate from top to bottom--------------------
        

        [snow_w(pS>0), snow_a(pS>0), water_flux] = infiltrateTop2Bottom(snow_i(pS>0), snow_w(pS>0), snow_a(pS>0), poreSpace(pS>0), maxWater, water_flux);
    
   
  
%------------infiltrate bottom to top----------------------

 % energyC3=sum(T.*c_temp.*(snow_i+snow_w+snow_a))+(sum(snow_w)+water_flux).*3.34e8;


        

        if water_flux>0 
            [snow_w(pS>0), snow_a(pS>0), runoff] = infiltrateBottom2Top(snow_i(pS>0), snow_w(pS>0), snow_a(pS>0), water_flux);
        end
    
  % energyC4=sum(T.*c_temp.*(snow_i+snow_w+snow_a))+(sum(snow_w)+runoff).*3.34e8;


%----------shift energy from water to T due to refreezing---------
    end
    
    [T, snow_i, snow_w] = refreeze(T, snow_i, snow_w, snow_a, c_temp, PARA);
       
  %  energyC5=sum(T.*c_temp.*(snow_i+snow_w+snow_a))+(sum(snow_w)+runoff).*3.34e8;
    
   % if abs((energyC5-energyC)./energyC)>0.01
    %    energyC
     %   energyC2
        
      %  energyC5
end
    
  
    
    


  
