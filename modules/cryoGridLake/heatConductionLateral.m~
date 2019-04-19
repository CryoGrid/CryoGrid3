function [dE_dt_cond dE_dt_lateral T_lateral]=heatConductionLateral(dE_dt_cond,T,k_temp,GRID,REF,FLAKE,t)

[~, idx]=min(abs(REF.load.OUT.timestamp-t));
T_lateral=REF.load.OUT.cryoGrid3(:,idx);

%Estimated lateral heat flux. This procedure assumes that the thermal regime
%of the ground is not affected by the lake in a distance of about 2 twice 
%the diameter of the lake. Assuming a circular lake shape, the lake has 
%an effetive cross section with the soil domain of pi. Thus, the lateral
%heat flux per m^2 scales with (3.14.*FLAKE.fetch.*GRID.general.K_delta)./(3.14./4.*FLAKE.fetch.^2)

dE_dt_lateral = k_temp.*(T_lateral-T)./(2.*FLAKE.fetch) .* 4.*GRID.general.K_delta./FLAKE.fetch;
dE_dt_lateral(~GRID.soil.cT_domain | T>0)=0; %exclude talik from lateral heat flux

dE_dt_cond=dE_dt_cond+dE_dt_lateral;

