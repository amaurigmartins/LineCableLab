function [Zg_self]=Z_der_slf(h,m_g,sigma_g,omega,con)
% Function for the Self Earth Impedance by Deri's Approximate
% Formula

% Inputs
% h1      : height of conductor [m]
% e_g     : permittivity of earth (originally not used) [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Z_self: Self earth impedance [Ohm/m]


m0=4*pi*1e-7;
p=1./sqrt(1i*omega.*m_g*sigma_g);
%p=1./sqrt(j*omega.*m_g.*(sigma_g+j*omega.*e_g));
Zg_self=zeros(con,con);
for k=1:1:con
    ter1_self=log((h(1,k)+p)/h(1,k));
    Zg_self(k,k)=1i*omega.*m0./(2*pi).*ter1_self;
end