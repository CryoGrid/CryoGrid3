<<<<<<< HEAD
function [T, snow_i, snow_w, snow_a]=melt(T, snow_i, snow_w, snow_a, poreSpace, c_temp)

L=3.34e8;

pot_SWE=double(T>0).*T.*c_temp.*(snow_i+snow_w+snow_a)./L;


T=double(T<=0).*T;
delta_SWE=pot_SWE;

if sum(pot_SWE>snow_i)~=0% in one cell more energy  than needed to melt entire cell
    for i=1:size(T,1)-1
        delta_SWE(i)=min([snow_i(i) pot_SWE(i)]);
        SWEres=pot_SWE(i) - delta_SWE(i);
        if (snow_i(i+1)+snow_w(i+1)+snow_a(i+1))~=0
            T(i+1)=T(i+1) + SWEres.*L./(c_temp(i+1).*(snow_i(i+1)+snow_w(i+1)+snow_a(i+1)));
        end
        pot_SWE(i+1)=pot_SWE(i+1) + double(T(i+1)>0).*T(i+1).*c_temp(i+1).*(snow_i(i+1)+snow_w(i+1)+snow_a(i+1))./L;   
        T(i+1)=double(T(i+1)<=0).*T(i+1);
    end
    i=size(T,1);
    delta_SWE(i)=min([snow_i(i) pot_SWE(i)]);
end
%T(find(isnan(T(:,1))==1),1)=0;

    
%delta_SWE
% delta_SWE = min([snow_i  double(T>0).*T.*c_temp.*(snow_i+snow_w+snow_a)./L]');   %Energy conserving since sensible heat term is only excess energy from SEB+conduction
% delta_SWE=delta_SWE';


snow_i=snow_i - delta_SWE; %melting only changes SWE, not density theta_s

snow_w=snow_w + delta_SWE;
%snow_a=(poreSpace.*snow_i + poreSpace.*snow_w - snow_w)./(1-poreSpace);  %pore space stays stays constant
=======
function [T, snow_i, snow_w, snow_a]=melt(T, snow_i, snow_w, snow_a, poreSpace, c_temp, PARA)

    L=PARA.constants.L_sl.*PARA.constants.rho_w; %3.34e8;

	pot_SWE=double(T>0).*T.*c_temp.*(snow_i+snow_w+snow_a)./L;


	T=double(T<=0).*T;
	delta_SWE=pot_SWE;

	if sum(pot_SWE>snow_i)~=0% in one cell more energy  than needed to melt entire cell
		for i=1:size(T,1)-1
		    delta_SWE(i)=min([snow_i(i) pot_SWE(i)]);
		    SWEres=pot_SWE(i) - delta_SWE(i);
		    if (snow_i(i+1)+snow_w(i+1)+snow_a(i+1))~=0
		        T(i+1)=T(i+1) + SWEres.*L./(c_temp(i+1).*(snow_i(i+1)+snow_w(i+1)+snow_a(i+1)));
		    end
		    pot_SWE(i+1)=pot_SWE(i+1) + double(T(i+1)>0).*T(i+1).*c_temp(i+1).*(snow_i(i+1)+snow_w(i+1)+snow_a(i+1))./L;   
		    T(i+1)=double(T(i+1)<=0).*T(i+1);
		end
		i=size(T,1);
		delta_SWE(i)=min([snow_i(i) pot_SWE(i)]);
	end
	%T(find(isnan(T(:,1))==1),1)=0;

		
	%delta_SWE
	% delta_SWE = min([snow_i  double(T>0).*T.*c_temp.*(snow_i+snow_w+snow_a)./L]');   %Energy conserving since sensible heat term is only excess energy from SEB+conduction
	% delta_SWE=delta_SWE';


	snow_i=snow_i - delta_SWE; %melting only changes SWE, not density theta_s

	snow_w=snow_w + delta_SWE;
	%snow_a=(poreSpace.*snow_i + poreSpace.*snow_w - snow_w)./(1-poreSpace);  %pore space stays stays constant
>>>>>>> origin/xice_mpi_polygon_TC

