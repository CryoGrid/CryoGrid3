function steadyState=steadyState(Tsurf,  Q, conductivity, cT_grid, K_grid,  K_frozen, K_thawed, arraySizeT)

posZero=find(cT_grid(:,1)>0);
posZero=posZero(1,1);

steadyState=zeros(size(cT_grid,1),1)-1;

steadyState(posZero)=Tsurf;
for i=posZero+1:size(steadyState,1)
    a=(steadyState(i-1,1)-K_frozen(i,1))./(K_thawed(i,1)-K_frozen(i,1))*(arraySizeT-2)+1; %T and c information live on same grid
posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));

   
    steadyState(i,1)=steadyState(i-1,1)+ 1/conductivity(i,posT)*Q*(cT_grid(i,1)-cT_grid(i-1,1));
end
