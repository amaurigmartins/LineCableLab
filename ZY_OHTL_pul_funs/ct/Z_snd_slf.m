% Z_snd_slf.m
function [Zg_self]=Z_snd_slf(num,h,e_g,m_g,sigma_g,omega)
% This function Calculates the Self Earth Impedance by Sunde's
% Approxiamtion

% input: num,h,e_g,_m_g,sigma_g,omega
% output: Z_self

Zg_self=zeros(num,num);
m0=4*pi*1e-7;
gamma_g=sqrt(1i*omega.*m_g.*(sigma_g+1i*omega.*e_g));

for k=1:1:num
   % Self Impedance
   ar_self=1+gamma_g.*h(1,k);
   par_self=gamma_g.*h(1,k);

    Zg_self(k,k)=1i*omega.*m0/(2*pi).*log(ar_self./par_self);
end

%save A.mat ar_self par_self