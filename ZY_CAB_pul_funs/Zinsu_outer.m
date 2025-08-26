% Zinsu_outer.m
function [zinsu]=Zinsu_outer(radius_ex,radius_in,mu_r,omega)
% Calculation of impedance of the insulation layer

% Constant Variables
m0=4*pi*1e-7;
mu_ins=m0*mu_r;


% Approximated skin effect
if (~radius_in)
    radius_in=eps;
end

zinsu = 1i.*omega.*mu_ins.*log(radius_ex./radius_in)/(2*pi);

end
