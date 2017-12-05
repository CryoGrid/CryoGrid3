function u_ini = initialTemperature(Temp_ini,z,z_help,capacity,conductivity,T_f,T_th,array_size)
Q_geo=-53e-3;
u_ini=interp1(Temp_ini(:,1),Temp_ini(:,2),z,'linear');
last_interp=find(~isnan(u_ini),1,'last');
dz=diff(z);
for i=last_interp:1:size(z,1)-1
   [c, k_help_z]=soil_thermalprop(u_ini([i i]),z(i:i+1),z_help(i:i+2),capacity(i:i+1,:),conductivity(i:i+1,:),T_f,T_th,array_size);
    u_ini(i+1)=-Q_geo.*dz(i)./k_help_z(2) + u_ini(i);
end
u_ini(isnan(u_ini))=u_ini(find(~isnan(u_ini),1,'first'));
