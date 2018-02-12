function DW = waterDensity(T)

%calulated following the CIPM formula
a1=-3.983035;
a2=301.797;
a3=522528.9;
a4=69.34881;
a5=999.974950;

DW=a5*(1- ((T+a1).^2.*(T+a2))./(a3.*(T+a4))); %[kg/m³]