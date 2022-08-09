function [sigma,rho,er] = soilFD_VisAl_fun(rho_100,f)

% Visacro & Alipio


sigma_eff_100=1/rho_100;

er=1.3+7.6*10^3.*f.^(-0.4);
sigma=abs(sigma_eff_100+1.2*10^(-6).*sigma_eff_100^0.27.*((f-100)).^0.65);
rho=1./sigma;





