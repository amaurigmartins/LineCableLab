function [sigma,rho,er] = soilFD_VisPor_fun(rho_100,f)

% Visacro & Portela


sigma_eff_100=1/rho_100;

er=2.34*10^6*sigma_eff_100^0.535.*f.^(-0.597);
sigma=sigma_eff_100.*(f./100).^0.072;
rho=1./sigma;





