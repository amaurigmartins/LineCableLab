function [Zg_self]=Z_alvadetan_slf(h,e_g,m_g,sigma_g,omega,con)
% Function for the Self Earth Impedance by Alvarado - Betancourt Method
% Approximate Formula

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

Zg_self=zeros(con,con);

for k=1:1:con

    Js=0.5*log((h(1,k)+p)./h(1,k))-1./24.*(1./(1+h(1,k)./p).^3);

    Zg_self(k,k)=1i*omega.*m0./(pi).*Js;
end