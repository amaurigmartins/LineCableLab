% Z_pet_mut.m
function [Zg_mutual]=Z_pet_mut(h1,h2,d,e_g,m_g,sigma_g,omega);
% This function Calculates the Mutual Earth Impedance by Semlyen's
% Approxiamtion

% input: h1,h2,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
e0=8.854187817e-12;
b=sqrt(j*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));

% Mutual Impedance
dp=sqrt((h1+h2+2./b).^2+d.^2);
d2=sqrt(d.^2+(h1+h2).^2);
term=log(dp./d2);


Zg_mutual=j*omega.*m0./(2*pi).*term;
