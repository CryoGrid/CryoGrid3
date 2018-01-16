function [Lambda]=conductionDevries(w_c,i_c,min_c,org_c,T,Xo,qo)
% function calculates soil heat conductivity using volumetric contents of water, ice,
% organic and minerals.
% Based on de Vries (1963) / Campbell et al. (1994)
% Input parameters: w_c:=water content
%                   i_c:=ice content
%                   min_c:=mineral content
%                   org_c:=organic content
%                   T:=soil temperature [°C]
%                   Xo:=cutoff water content for liquid recirculation [m³/m³]
%                   qo:=power for liquid recirculation function
%--------------------------------------------------------------------------


warning off all
P=1-min_c-org_c;

L(1,:)=0.025;       %air [Hillel(1982)]
L(2,:)=0.57;        %water [Hillel(1982)]
L(3,:)=2.2;         %ice [Hillel(1982)]
L(4,:)=0.25;        %organic [Hillel(1982)]
L(5,:)=3.8;%7.7;%y2.9;         %mineral [Hillel(1982)]

%calaculated conductivties [Campbell et al. (1994)]
%L(1,:)=0.024+7.73*10^-5*T-2.6*10^-8*T^2; %air
%L(2,:)=0.554+2.24*10^-3*T-9.87*10^-6*T^2;%water

%--------------------------------------------------------------------------
%conductivity function due to variation in water-air content (no ice)
wc=[0:0.01:P];
ac=P-wc;
m=size(wc,2);
space=zeros(m*3, 3);

for i=1:size(wc,2)
    X(1,:)=ac(i);        %air
    X(2,:)=wc(i);        %water
    X(3,:)=0.0;          %ice
    X(4,:)=org_c;        %organic
    X(5,:)=min_c;        %mineral
    %----------------------------------------------------------------------
    [kwa(i)]=conductDV63_aw(X,L,T,Xo,qo,[]);
    KK(i,1)=wc(i);
    KK(i,2)=0;
    KK(i,3)=kwa(i);
    space(i,1)=wc(i);
    space(i,3)=kwa(i);
end
%plot(space(1:m,1), space(1:m,3))
KK(end,:)=[];
l=size(KK,1);
%--------------------------------------------------------------------------
%conductivity function due to variation in ice-air content (no watr)
ic=[0:0.01:P];
ac=P-ic;

for i=1:size(ic,2)
    X(1,:)=ac(i);        %air
    X(2,:)=0.0;          %water
    X(3,:)=ic(i);        %ice  
    X(4,:)=org_c;        %organic
    X(5,:)=min_c;        %mineral
    %----------------------------------------------------------------------
    [kia(i)]=conductDV63_ia(X,L,T,Xo,qo,[]);
    KK(l+i,1)=0;
    KK(l+i,2)=ic(i);
    KK(l+i,3)=kia(i);
    space(m+i,2)=ic(i);
    space(m+i,3)=kia(i);
end

KK(end,:)=[];
l=size(KK,1);
%-------------------------------------------------------------------------
%conductivity function due to variation in ice-water content (no air)
ic=[0:0.01:P];
wc=P-ic;

for i=1:size(wc,2)
    X(1,:)=0.0;          %air
    X(2,:)=wc(i);        %water
    X(3,:)=ic(i);        %ice  
    X(4,:)=org_c;        %organic
    X(5,:)=min_c;        %mineral
    %----------------------------------------------------------------------
    [kiw(i)]=conductDV63_iw(X,L,T,Xo,qo,[]);
    KK(l+i,1)=wc(i);
    KK(l+i,2)=ic(i);
    KK(l+i,3)=kiw(i);
    space(2*m+i,1)=wc(i);
    space(2*m+i,2)=ic(i);
    space(2*m+i,3)=kiw(i);

end

KK(end,:)=[0 P kiw(end)];
KK(1,:)=[];
%--------------------------------------------------------------------------
%interpolate boundary functions and calculate searched conductivity

%Lamda = griddata(KK(:,1),KK(:,2),KK(:,3),XI,YI,'v4');

max_size1=max([size(1-min_c-org_c,1) size(w_c,1) size(i_c,1)]);
max_size2=max([size(1-min_c-org_c,2) size(w_c,2) size(i_c,2)]);
if max_size1>max_size2 %column vector
    a_c=zeros(max_size1,1)+1-min_c-org_c-w_c-i_c;
    w_c=zeros(max_size1,1)+w_c;
    i_c=zeros(max_size1,1)+i_c;
else %row vector
    a_c=zeros(1,max_size2)+1-min_c-org_c-w_c-i_c;
    w_c=zeros(1, max_size2)+w_c;
    i_c=zeros(1, max_size2)+i_c;
end
Lambda=griddata(space(:,1), space(:,2), space(:,3),w_c, i_c, 'v4');
%F=TriScatteredInterp(space(:,1), space(:,2), space(:,3), space(:,4));
%Lambda=F(a_c, w_c, i_c)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [k]=conductDV63_aw(X,L,T,Xo,qo,fw)
%function calculates heat conductivity based on de Vries (1963)
%Input parameter:   X:=fraction of soil material
%                   L:=individual conductivities of soil components
%                   Xo:=cutoff water content for liquid recirculation [m³/m³]
%                   qo:=power for liquid recirculation function
%                   fw:=weighting function (if empty else fw=1 for saturated
%                   soil or fw=0 for dry soils 
%--------------------------------------------------------------------------
%g:=shape factores of soil components (set to spherical gx=gy=gz)
g(1,:)=1/3;
g(2,:)=1/3;
g(3,:)=1/3;
%power for liquid recirculation function
qo=qo;
%cutoff water content for liquid recirculation [m³/m³]
Xo=Xo;
%calculate Ly:=conductivity of matrix material 
q=qo*((T+273.15)/303)^2;
if isempty(fw)==1
    fw=1/(1+(X(2,:)/Xo)^-q);
end
Ly=L(1,:)+fw*(L(2,:)-L(1,:));
%calculate k:=effectiv conductivity 
for i=1:5
    for j=1:3
        K1(i,j)=1/(1+(L(i,1)/Ly-1)*g(j,1));
    end
end
K=1/3*sum(K1')';
k=sum(K.*L.*X)/sum(K.*X);


function [k]=conductDV63_ia(X,L,T,Xo,qo,fw)
%function calculates heat conductivity based on de Vries (1963)
%Input parameter:   X:=fraction of soil material
%                   L:=individual conductivities of soil components
%                   Xo:=cutoff water content for liquid recirculation [m³/m³]
%                   qo:=power for liquid recirculation function
%                   fw:=weighting function (if empty else fw=1 for saturated
%                   soil or fw=0 for dry soils 
%--------------------------------------------------------------------------
%g:=shape factores of soil components (set to spherical gx=gy=gz)
g(1,:)=1/3;
g(2,:)=1/3;
g(3,:)=1/3;
%power for liquid recirculation function
qo=qo;
%cutoff water content for liquid recirculation [m³/m³]
Xo=Xo;
%calculate Ly:=conductivity of matrix material 
q=qo*((T+273.15)/303)^2;
if isempty(fw)==1
    fw=1/(1+(X(3,:)/Xo)^-q);
end
Ly=L(1,:)+fw*(L(3,:)-L(1,:));
%calculate k:=effectiv conductivity 
for i=1:5
    for j=1:3
        K1(i,j)=1/(1+(L(i,1)/Ly-1)*g(j,1));
    end
end
K=1/3*sum(K1')';
k=sum(K.*L.*X)/sum(K.*X);


function [k]=conductDV63_iw(X,L,T,Xo,qo,fw)
%function calculates heat conductivity based on de Vries (1963)
%Input parameter:   X:=fraction of soil material
%                   L:=individual conductivities of soil components
%                   Xo:=cutoff water content for liquid recirculation [m³/m³]
%                   qo:=power for liquid recirculation function
%                   fw:=weighting function (if empty else fw=1 for saturated
%                   soil or fw=0 for dry soils 
%--------------------------------------------------------------------------
%g:=shape factores of soil components (set to spherical gx=gy=gz)
g(1,:)=1/3;
g(2,:)=1/3;
g(3,:)=1/3;
%power for liquid recirculation function
qo=qo;
%cutoff water content for liquid recirculation [m³/m³]
Xo=Xo;
%calculate Ly:=conductivity of matrix material 
q=qo*((T+273.15)/303)^2;
if isempty(fw)==1
    fw=1/(1+(X(2,:)/Xo)^-q);
end
Ly=L(3,:)+fw*(L(2,:)-L(3,:));
%calculate k:=effectiv conductivity 
for i=1:5
    for j=1:3
        K1(i,j)=1/(1+(L(i,1)/Ly-1)*g(j,1));
    end
end
K=1/3*sum(K1')';
k=sum(K.*L.*X)/sum(K.*X);
