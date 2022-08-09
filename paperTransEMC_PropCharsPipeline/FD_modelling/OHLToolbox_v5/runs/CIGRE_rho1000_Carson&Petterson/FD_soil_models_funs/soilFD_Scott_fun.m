function [sigma,rho,er] = soilFD_Scott_fun(rho_100,f)

% Scott

sigma_100=1000/rho_100;


D1=5.491+0.946.*log10(sigma_100)-1.097*log10(f)+0.069.*(log10(sigma_100)).^2-0.114.*log10(sigma_100).*log10(f)+0.067.*(log10(f)).^2;
D2=0.028+1.098.*log10(sigma_100)-0.068*log10(f)+0.036.*(log10(sigma_100)).^2-0.046.*log10(sigma_100).*log10(f)+0.018.*(log10(f)).^2;

er=10.^D1;
sigma=10.^D2;
rho=1000./sigma;

