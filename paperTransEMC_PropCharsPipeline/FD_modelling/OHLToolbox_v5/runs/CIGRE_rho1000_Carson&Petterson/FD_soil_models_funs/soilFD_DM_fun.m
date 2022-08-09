function [sigma,rho,er] = soilFD_DM_fun(rho_LF,er_HF_dry,f)

% Datsios - Mikropoulos
sigma_LF=1E4/rho_LF;
f_LF=42;
% er calculation
K=0.537*sigma_LF^0.16;
er_3k_dry=2.9*er_HF_dry-3.8;
er_HF=(1.24*sigma_LF^0.415)*er_HF_dry;
er_3k=(4.00*sigma_LF^0.463)*er_3k_dry;
er=er_HF+(3000./f).^K*(er_3k-er_HF);
%er=er_HF+(3000./f).^K*(er_3k-K);
indxlow3k=f(f(:,1)<3000);
[I,J] = find(f==3000);
er(1:length(indxlow3k),1)=er(I);

% sigma calculation
sigma_HF=sigma_LF*(1+0.65/(sigma_LF)^0.57);

sigma=sigma_HF.*f.*1E-6+(f_LF-f_LF*(f-f_LF).*1E-6).*(sigma_LF/f_LF-sigma_HF*1E-6); % (uS/cm)
rho=1E4./sigma; %(Ohmm);

% figure (1)
% subplot(2,1,1)
% loglog(f,sigma)
% xlabel('frequency (Hz)')
% ylabel('\muS/cm')
% subplot(2,1,2)
% semilogx(f,rho)
% xlabel('frequency (Hz)')
% ylabel('Ohm/m')
% 
% figure (2)
% loglog(f,er)
% xlabel('frequency (Hz)')
% ylabel('er')
