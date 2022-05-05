function [sigmaeff,rhoeff,er] = soilFD_AlVis_fun(rho_100,f)

% Alipio Visacro

e0=8.85418782E-12;
sigma_100=1/rho_100;
sigma_eff_100=1000/rho_100;
gamma1=0.54;
er_inf=12;
h1=1.26*sigma_eff_100^(-0.73);

er=er_inf+(tan(pi*gamma1/2)*1E-3)/(2*pi*e0*10^(6*gamma1))*sigma_eff_100*h1.*f.^(gamma1-1);
sigmaeff=sigma_eff_100+sigma_eff_100*h1.*(f./1E6).^gamma1;

rhoeff=1000./sigmaeff;




