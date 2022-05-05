function [sigma,rho,er] = soilFD_LongSmith_fun(rho_100,er_HF,f)

% Longmire & Smith
an=load('coeff.txt');
e0=8.85418782E-12;
sigma_100=1/rho_100;
siz=length(f);

er_sum=zeros(siz,1);
sigma_sum=zeros(siz,1);

for n=1:1:13
    [e_term,s_term] = sum_fun(sigma_100,f,an(n),n);
    er_sum=er_sum+e_term;
    sigma_sum=sigma_sum+s_term;
end


er=er_HF+er_sum;
%er=er_HF+sigma_sum;
sigma=sigma_100+2*pi.*f*e0.*sigma_sum;
rho=1./sigma;


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
