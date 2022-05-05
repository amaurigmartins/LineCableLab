% Z_pet_slf.m
function [Yg_self,log_ins]=Y_pet_slf_imperf(h1,e_g,m_g,sigma_g,omega)
% This function Calculates the self potential coefficient (imperfect ground - Pettersson's Approximation)

% inputs:
% h1: height of conductor
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Pg_self

m0=4*pi*1e-7;
e0=8.854187817e-12;
b=sqrt(1i*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));

n2=sqrt((1i.*omega*e_g+sigma_g)./(1i.*omega.*e0));
n1=1;
n=n2./n1;

% Self potential coefficient 
log_term_ar=n.^2+1;
log_term_par=2.*b.*h1;
log_term_ins=log_term_ar./log_term_par;
log_ins=(1+log_term_ins);
log_term=log(log_ins);

Yg_self=2./(n.^2+1).*log_term;
