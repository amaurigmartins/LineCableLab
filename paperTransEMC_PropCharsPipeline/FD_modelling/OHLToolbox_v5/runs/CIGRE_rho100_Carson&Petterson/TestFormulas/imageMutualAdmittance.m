function Ym_image = imageMutualAdmittance(f,hi,hj,x)

w = 2*pi*f;

% Constants
e0 = 8.854187817e-12; % Vacuum permittivity
% Geometry formulas
H=hi+hj;
dij = sqrt((hi-hj)^2+x^2);
Dij = sqrt(H^2+x^2);


Pkm = (1/(2*pi*e0))*log(Dij/dij);
Ym_image = 1i*w*1/Pkm;

end