function res= psi_H(zeta1, zeta2)

if zeta1<=0
    
    res= 1.9.*atanh((1 - 11.6.*zeta1).^0.5) + log(zeta1) - (1.9.*atanh((1 - 11.6.*zeta2).^0.5) + log(zeta2));
   
    
   % res2=quad(@(z) (1-0.95.*(1-11.6.*z).^(-0.5))./z, zeta2, zeta1); 
else
    res=((-5 + 5^0.5).*log(-3 + 5^0.5- 2.*zeta1) - (5 + 5^0.5).*log(3 + 5^0.5 + 2.*zeta1))/2  - (((-5 + 5^0.5).*log(-3 + 5^0.5- 2.*zeta2) - (5 + 5^0.5).*log(3 + 5^0.5 + 2.*zeta2))/2);
    
  %  res2=quad(@(z) (1-(1+(5.*z.*(1+z))./(1+3.*z+z.^2)))./z, zeta2, zeta1); 
end


%res=quad(@(z) (1-(0.95+7.8*z))./z, zeta2, zeta1); 