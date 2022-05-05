% Z_snd_slf.m
function [Zg_self]=Z_snd_slf(h,e_g,m_g,sigma_g,omega,con)
% Function for the Self Earth Impedance by Pettersson's Approximate
% Formula

% Inputs
% h1      : height of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Zg_self : Self earth impedance [Ohm/m]

Zg_self=zeros(con,con);

m0=4*pi*1e-7;
gamma_g=sqrt(1i*omega.*m_g.*(sigma_g+1i*omega.*e_g));

for k=1:1:con
    % Self Impedance
    ar_self=1+gamma_g.*h(1,k);
    par_self=gamma_g.*h(1,k);

    Zg_self(k,k)=1i*omega.*m0/(2*pi).*log(ar_self./par_self);
end