function Ys_Wise = wiseSelfAdmittance(f,h,ri,rho,mu_rg,e_rg)
% This function Calculates the Self Generalized Earth Admittance by Wise's
%
% INPUTS:
%       f: frequency [Hz]
%       h: height of conductor [m]
%       ri: outer conductor radius [m]
%       rho: earth resistivity [ohm.m]
%       e_rg: relative permittivity of earth [p.u.]
%       mu_rg: relative permeability of earth [p.u.]
%
% OUTPUT: Ys_Wise : self admittance [S/m]


w = 2*pi*f;

% Constants
e0 = 8.854187817e-12; % Vacuum permittivity
mu0 = 4*pi*1e-7; % Vacuum permeability
k0 = w*(sqrt(mu0*e0));
gama0 = 1i*w*sqrt(mu0*e0);

% Geometry formulas
H=2*h;
dij = ri;
Dij = sqrt(H^2);

% Ground parameters
sigma_g = 1/rho;
eg = e0*e_rg;
mu_g = mu0*mu_rg;
gama_g = sqrt(1i*w*mu0*(sigma_g + 1i*w*eg));

% Adimittance formula:
y_lambda =@(lambda) (((2*mu_g*(gama0^2)*(mu0*lambda + (sqrt(lambda^2 + gama_g^2 + k0^2))*mu_g)*exp(-H*lambda))/(((sqrt(lambda^2 + gama_g^2 + k0^2))*mu0 + lambda*mu_g)*((sqrt(lambda^2 + gama_g^2 + k0^2))*(gama0^2)*mu_g + lambda*(gama_g^2)*mu0))));
Qm = integral(y_lambda,0,Inf,'ArrayValued',true);
Peij = 1/(2*pi*e0)*(log(Dij/dij)+Qm);

Ys_Wise = 1i*w*(Peij^(-1));

end
