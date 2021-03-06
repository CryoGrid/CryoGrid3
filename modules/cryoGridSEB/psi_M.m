function res= psi_M(zeta1, zeta2)

if zeta1<=0
    
    res=-2.*atan((1 - 19.*zeta1).^(1/4)) + 2.* log(1 + (1 - 19.*zeta1).^(1/4)) + log(1 + (1 - 19.*zeta1).^0.5) - ...
        (-2.*atan((1 - 19.*zeta2).^(1/4)) + 2.* log(1 + (1 - 19.*zeta2).^(1/4)) + log(1 + (1 - 19.*zeta2).^0.5));
    
   % res2=quad(@(z) (1-(1-19.*z).^(-0.25))./z, zeta2, zeta1); 
    
else
    
    res=-19.5.*(1 + zeta1).^(1/3) - 7.5367.*atan(0.57735 - 1.72489.*(1 + zeta1).^(1/3)) + 4.35131.*log(3+4.4814.*(1+zeta1).^(1/3)) - 2.17566.*log(3 - 4.4814.*(1 + zeta1).^(1/3) + 6.69433.*(1 + zeta1).^(2/3)) - ...
        (-19.5.*(1 + zeta2).^(1/3) - 7.5367.*atan(0.57735 - 1.72489.*(1 + zeta2).^(1/3)) + 4.35131.*log(3+4.4814.*(1+zeta2).^(1/3)) - 2.17566.*log(3 - 4.4814.*(1 + zeta2).^(1/3) + 6.69433.*(1 + zeta2).^(2/3))) ;
    
   % res2=quad(@(z) (1-(1+(6.5.*z.*(1+z).^(1/3))./(1.3+z)))./z, zeta2, zeta1); 
end

 %res=quad(@(z) (1-(1+6*z))./z, zeta2, zeta1); 