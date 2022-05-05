% Z_snd_mut.m
function [Zg_mutual]=Z_snd2_mut(h1,h2,d,e_g,m_g,sigma_g,omega);
% This function Calculates the Mutual Earth Impedance by Sunde's
% Approxiamtion

% input: h1,h2,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
e0=8.854187817e-12;

gamma_g=j*omega.*m_g.*(sigma_g+j*omega.*e_g);
gamma_a=j*omega.*m0.*(j*omega.*e0);
gamma_eq=sqrt(gamma_g-gamma_a);

% Mutual Impedance
ar_mut=(1+gamma_eq.*((h1+h2)/2)).^2+(gamma_eq.*(d/2)).^2;
par_mut=(gamma_eq.*((h1+h2)/2)).^2+(gamma_eq.*(d/2)).^2;

Zg_mutual=j*omega.*m0/(4*pi).*log(ar_mut./par_mut);