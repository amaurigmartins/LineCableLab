function [sigma,rho,er] = soilFD_Cigre_fun(rho_100,f)

% Visacro & Alipio


sigma_eff_100=1/rho_100;

er=12+9.5*10^4.*sigma_eff_100^0.27.*f.^(-0.46);
sigma=sigma_eff_100+4.7*10^(-6).*sigma_eff_100^0.27.*(f).^0.54;
rho=1./sigma;





