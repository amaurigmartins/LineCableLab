% Potcoef_outer.m
function [pm]=Potcoef_outer(radius_ex,radius_in,eps_r,omega,rho_ins)
% Calculation of the potential coefficient for a tubular insulation layer

% Parse optional input
if nargin < 5
    rho_ins = Inf;
end

% Constant Variables
e0=8.854187817e-12;
eps_ins=e0.*eps_r + 1./(1i.*omega.*rho_ins);

pm = log(radius_ex./radius_in)./(2.*pi.*eps_ins);

end
