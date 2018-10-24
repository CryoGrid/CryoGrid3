function [c_temp, k_eff, lwc_temp]=readThermalParameters(T, GRID, PARA)


cT_frozen = GRID.soil.cT_frozen;
cT_thawed = GRID.soil.cT_thawed;
capacity = GRID.soil.capacity;
K_frozen = GRID.soil.K_frozen; 
K_thawed = GRID.soil.K_thawed;
conductivity = GRID.soil.conductivity; 
arraySizeT = PARA.technical.arraySizeT;
liquidWaterContent = GRID.soil.liquidWaterContent;
                                               

a=(T-cT_frozen)./(cT_thawed-cT_frozen)*(arraySizeT-2)+1; %T and c information live on same grid
posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));
posT(posT==0)=1;
posT(isnan(posT))=arraySizeT;


r=[1:size(capacity,1)]';
c=posT;
indices=(size(capacity,1))*(c-1)+r;
c_temp=capacity(indices);
lwc_temp=liquidWaterContent(indices);

a=(T-K_frozen)./(K_thawed-K_frozen)*(arraySizeT-2)+1;
posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));
posT(posT==0)=1;
posT(isnan(posT))=arraySizeT;


r=[1:size(conductivity,1)]';
c=posT;
indices=(size(conductivity,1))*(c-1)+r;
k_eff=conductivity(indices);
