% Z_pet_slf.m
function [Zg_pet_self]=Z_pet_slf(num,h,e_g,m_g,sigma_g,omega)
% This function Calculates the Self Earth Impedance by Deri's
% Approxiamtion

% input: num,h,e_g,_m_g,sigma_g,omega
% output: Zg_pet_self
m0=4*pi*1e-7;
e0=8.854187817e-12;
b=sqrt(1i*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));

Zg_pet_self=zeros(num,num);

for k=1:1:num
    
  % Self Impedance
  ar_self=1+b*h(1,k);
  par_self=b*h(1,k);

  Zg_pet_self(k,k)=j*omega.*m0./(2*pi).*log(ar_self./par_self);

end



