% Z_snd_slf.m
function [Zg_self]=Z_snd2_slf(h1,e_g,m_g,sigma_g,omega);
% This function Calculates the Self Earth Impedance by Sunde's
% Approxiamtion

% input: h1,e_g,_m_g,sigma_g,omega
% output: Z_self
m0=4*pi*1e-7;
e0=8.854187817e-12;

gamma_g=j*omega.*m_g.*(sigma_g+j*omega.*e_g);
gamma_a=j*omega.*m0.*(j*omega.*e0);
gamma_eq=sqrt(gamma_g-gamma_a);

% Self Impedance
ar_self=1+gamma_eq.*h1;
par_self=gamma_eq.*h1;

Zg_self=j*omega.*m0/(2*pi).*log(ar_self./par_self);
