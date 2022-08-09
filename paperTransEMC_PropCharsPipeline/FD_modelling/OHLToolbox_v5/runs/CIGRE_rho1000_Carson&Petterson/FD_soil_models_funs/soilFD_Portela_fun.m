function [sigmaeff,rhoeff,er] = soilFD_Portela_fun(rho_100,f)

% Portela
e0=8.85418782E-12;
beta=0.1;
ap=0.72;
sigmaeff_100=1/rho_100;

sigmaeff=sigmaeff_100*1E6+beta.*(2*pi.*f).^ap;
sigmaeff=sigmaeff/100;
rhoeff=1E4./sigmaeff;

er=((beta.*tan(pi./2*ap)*1E-6)/e0).*(2*pi.*f).^(ap-1);


