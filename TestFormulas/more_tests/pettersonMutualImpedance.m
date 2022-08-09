function Zm_Petterson = pettersonMutualImpedance(f,hi,hj,x,rho,e_rg)
% This function Calculates the Pettersson’s Closed-Form Expressions for
% Impedance
%
% INPUTS:
%       f: frequency [Hz]
%       hi: height of conductor i [m]
%       hj: height of conductor j [m]
%       x: horizontal distance between conductors i and j [m]
%       rho: earth resistivity [ohm.m]
%       e_rg: relative permittivity of earth [p.u.]
%       mu_rg: relative permeability of earth [p.u.]
%
% OUTPUT: Zm_Petterson : mutual impedance between conductors i and j [ohm/m]


w = 2*pi*f;

% Constants
e0 = 8.854187817e-12; % Vacuum permittivity
mu0 = 4*pi*1e-7; % Vacuum permeability
gama0 = 1i*w*sqrt(mu0*e0);

% Geometry formulas
H=hi+hj;
dij = sqrt((hi-hj)^2+x^2);
Dij = sqrt(H^2+x^2);

% Ground parameters
sigma_g = 1/rho;
eg = e0*e_rg;
gama_g = sqrt(1i*w*mu0*(sigma_g + 1i*w*eg));

% Petterson's formula:

beta = sqrt(gama_g^2 - gama0^2);

Me = log((sqrt((H + (2/(beta))^2 + x^2)))/(Dij));

Zm_Petterson = 1i*w*mu0/(2*pi)*(log(Dij/dij)+Me);

end
