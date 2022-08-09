% Z_pet_slf.m
function [Zg_self]=Z_pet_slf(h1,e_g,m_g,sigma_g,omega);
% This function Calculates the Self Earth Impedance by Deri's
% Approxiamtion

% input: h1,e_g,_m_g,sigma_g,omega
% output: Z_self
m0=4*pi*1e-7;
e0=8.854187817e-12;
b=sqrt(j*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));


% Self Impedance
ar_self=1+b*h1;
par_self=b*h1;

Zg_self=j*omega.*m0./(2*pi).*log(ar_self./par_self);
