function [sigma,rho,er] = soilFD_M_fun(rho_DC,er_HF,f)

% Messier model

e0=8.85418782E-12;

sigma_DC=1/rho_DC;


er=er_HF+sqrt((sigma_DC.*er_HF)./(pi.*f.*e0));


sigma=sigma_DC+sqrt(4*pi.*f*sigma_DC*e0*er_HF);
rho=1./sigma;