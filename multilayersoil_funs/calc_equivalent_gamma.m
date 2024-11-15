function [gamma_eq, eps_eff, sigma_eff] = calc_equivalent_gamma(sigma, epsr, t, freq)
% calc_equivalent_gamma - Calculates the equivalent propagation constant
% for an N-layer earth structure using conductivity, relative permittivity,
% and layer thicknesses, using the EHEM formulation of Xue et al.
% (10.1109/TPWRD.2021.3049595). The effective conductivity extracted from
% the imaginary part of `gamma_eq` is numerically identical to the
% equivalent conductivity of Martins-Britto et al. (10.1109/TPWRD.2019.2930406).
%
% Arguments:
% sigma : vector of conductivities for each layer [S/m]
% epsr  : vector of relative permittivities for each layer
% t     : vector of thicknesses for each layer (N-1 elements, bottom layer is infinite)
% freq  : frequency of interest [Hz]
%
% Returns:
% gamma_eq : Equivalent propagation constant for the N-layer structure.
% eps_eff  : Effective permittivity
% sigma_eff: Effective conductivity

% Validate inputs
if length(sigma) ~= length(epsr)
    error('sigma and epsr must have the same length, representing the number of layers.');
end
if length(t) ~= length(sigma) - 1
    error('t must have N-1 elements, one for each inter-layer boundary.');
end
if size(sigma,2) > 1; sigma = sigma.'; end
if size(epsr,2) > 1; epsr = epsr.'; end
if size(t,2) > 1; t = t.'; end

mu0 = 4*pi*1e-7;
eps0= 8.854187817e-12;

omega = 2*pi*freq;
gamma0_sq = -mu0*eps0*omega.^2;
gamma1_sq = 1i*omega*mu0.*(1i*omega.*eps0.*epsr+sigma);

gamma_e=sqrt(gamma1_sq-gamma0_sq);


i=size(gamma_e,1);
while 1

    if i==size(gamma_e,1)

        num=(gamma_e(i-1,:)+gamma_e(i,:)-(gamma_e(i-1,:)-gamma_e(i,:)).*exp(-2*t(i-1)*gamma_e(i-1,:)));
        den=(gamma_e(i-1,:)+gamma_e(i,:)+(gamma_e(i-1,:)-gamma_e(i,:)).*exp(-2*t(i-1)*gamma_e(i-1,:)));
        gamma_eq=gamma_e(i-1,:).*(num./den);

        if i==2
            break
        end

    elseif i~=[size(gamma_e,1),2]
        num=(gamma_e(i-1,:)+gamma_eq-(gamma_e(i-1,:)-gamma_eq).*exp(-2*t(i-1)*gamma_e(i-1,:)));
        den=(gamma_e(i-1,:)+gamma_eq+(gamma_e(i-1,:)-gamma_eq).*exp(-2*t(i-1)*gamma_e(i-1,:)));
        gamma_eq=gamma_e(i-1,:).*(num./den);
    elseif i==2
        num=(gamma_e(i-1,:)+gamma_eq-(gamma_e(i-1,:)-gamma_eq).*exp(-2*t(i-1)*gamma_e(i-1,:)));
        den=(gamma_e(i-1,:)+gamma_eq+(gamma_e(i-1,:)-gamma_eq).*exp(-2*t(i-1)*gamma_e(i-1,:)));
        gamma_eq=gamma_e(i-1,:).*(num./den);
        break
    end
    i=i-1;

end

gamma_eq_sq=gamma_eq.^2;

% Effective permittivity and conductivity
eps_eff = -real(gamma_eq_sq) / (omega^2 * mu0);
sigma_eff = imag(gamma_eq_sq) / (omega * mu0);

end