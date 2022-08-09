% Z_snd_slf.m
function [Zg_der_self]=Z_der_slf(num,h,e_g,m_g,sigma_g,omega)
% This function Calculates the Self Earth Impedance by Deri's
% Approxiamtion

% input: num,h1,e_g,_m_g,sigma_g,omega
% output: Z_der_self
m0=4*pi*1e-7;
p=1./sqrt(1i*omega.*m_g*sigma_g);

Zg_der_self=zeros(num,num);
%p=1./sqrt(j*omega.*m_g.*(sigma_g+j*omega.*e_g));

for k=1:1:num
% Self Impedance
ter1_self=log((h(1,k)+p)/h(1,k));


Zg_der_self(k,k)=1i*omega.*m_g./(2*pi).*ter1_self;
end
