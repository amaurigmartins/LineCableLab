% Z_snd_slf.m
function [Zg_sln_self]=Z_sln_slf(num,h,e_g,m_g,sigma_g,omega)
% This function Calculates the Self Earth Impedance by Semlyen's
% Approxiamtion

% input: h1,e_g,_m_g,sigma_g,omega
% output: Z_self
m0=4*pi*1e-7;

%Set a matrix Zsln_self with zeros

Zg_sln_self=zeros(num,num);

for k=1:1:num
    
    % Self Impedance
    ar_self=1i*omega.*m_g;
    par_self=sigma_g+1i*omega.*e_g;

    Zg_sln_self(k,k)=1/(2*pi*h(1,k))*sqrt(ar_self./par_self);
    
end

