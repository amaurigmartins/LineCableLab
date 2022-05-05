% Z_pet_slf.m
function [Yg_mut,dq]=Y_pet_mut_imperf(h1,h2,d,e_g,m_g,sigma_g,omega)
% This function Calculates the self potential coefficient (imperfect ground - Pettersson's Approximation)

% inputs:
% h1: height of conductor
% h2: height of conductor
% d12: distance of conductors
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Pg_mutual

% Mutual Impedance
m0=4*pi*1e-7;
e0=8.854187817e-12;

b=sqrt(1i*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));
n2=sqrt((1i.*omega*e_g+sigma_g)./(1i.*omega.*e0));
n1=1;
n=n2./n1;

d2=sqrt(d.^2+(h1+h2).^2);

term_dq=(n.^2+1)./b;
dq=-sqrt((h1+h2+term_dq).^2+d.^2);

temp=imag(dq);

if temp>0
    dq=sqrt((h1+h2+term_dq).^2+d.^2);
end

% Self Impedance
mut_term2=log(dq./d2);
mut_term1=2./(n.^2+1);

Yg_mut=mut_term1.*mut_term2;