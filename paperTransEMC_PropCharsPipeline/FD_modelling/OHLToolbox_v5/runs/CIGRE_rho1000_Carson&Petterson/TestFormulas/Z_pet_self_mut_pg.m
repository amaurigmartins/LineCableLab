% Z_pet_mut.m
function [Zpg_mat]=Z_pet_mut(h1,h2,h3,d12,d23,d13,r);
% This function Calculates the Mutual Earth Impedance by Semlyen's
% Approxiamtion

% input: h1,h2,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
e0=8.854187817e-12;


% Self Impedance
term_self=log((2*h1)./r);
Zpg11self=j*omega.*m0./(2*pi).*term_self;

term_self=log((2*h2)./r);
Zpg22self=j*omega.*m0./(2*pi).*term_self;

term_self=log((2*h3)./r);
Zpg33self=j*omega.*m0./(2*pi).*term_self;

% Mutual Impedance
d1=sqrt(d12.^2+(h1-h2).^2);
d2=sqrt(d12.^2+(h1+h2).^2);
term_mut=log(d2./d1);
Zpg12_mut=j*omega.*m0./(2*pi).*term_mut;

d1=sqrt(d13.^2+(h1-h3).^2);
d2=sqrt(d13.^2+(h1+h3).^2);
term_mut=log(d2./d1);
Zpg13_mut=j*omega.*m0./(2*pi).*term_mut;

d1=sqrt(d23.^2+(h2-h3).^2);
d2=sqrt(d23.^2+(h2+h3).^2);
term_mut=log(d2./d1);
Zpg23_mut=j*omega.*m0./(2*pi).*term_mut;
