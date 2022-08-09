function Zs_Wise = wiseSelfImpedance(f,h,ri,rho,mu_rg,e_rg)
% This function Calculates the Self Generalized Earth Impedance by Wise's
%
% INPUTS:
%       f: frequency [Hz]
%       h: height of conductor [m]
%       ri: outer conductor radius [m]
%       rho: earth resistivity [ohm.m]
%       e_rg: relative permittivity of earth [p.u.]
%       mu_rg: relative permeability of earth [p.u.]
%
% OUTPUT: Zs_Wise : self impedance [ohm/m]

w = 2*pi*f;

% Constants
e0 = 8.854187817e-12; % Vacuum permittivity
mu0 = 4*pi*1e-7; % Vacuum permeability
k0 = w*(sqrt(mu0*e0));

% Geometry formulas
H=2*h;
dij = ri;
Dij = sqrt(H^2);

% Ground parameters
sigma_g = 1/rho;
eg = e0*e_rg;
mu_g = mu0*mu_rg;
gama_g = sqrt(1i*w*mu0*(sigma_g + 1i*w*eg));


z_lambda =@(lambda) ((2*mu_g*(exp(-H*lambda)))/((sqrt(lambda^2 + gama_g^2 + k0^2))*mu0 + lambda*mu_g));
Jm = integral(z_lambda,0,Inf,'ArrayValued',true);
Zs_Wise = 1i*w*mu0/(2*pi)*(log(Dij/dij)+Jm);

end