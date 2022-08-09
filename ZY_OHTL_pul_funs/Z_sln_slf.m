% Z_snd_slf.m
function [Zg_self]=Z_sln_slf(h,e_g,m_g,sigma_g,omega,con)
% Function for the Self Earth Impedance by Semlyen's Approximate
% Formula

% Inputs
% h1      : height of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Zg_self: Self earth impedance [Ohm/m]

Zg_self=zeros(con,con);

% Self Impedance
for k=1:1:con
    ar_self=1i*omega.*m_g;
    par_self=sigma_g+1i*omega.*e_g;

    Zg_self(k,k)=1/(2*pi*h(1,k))*sqrt(ar_self./par_self);
end
