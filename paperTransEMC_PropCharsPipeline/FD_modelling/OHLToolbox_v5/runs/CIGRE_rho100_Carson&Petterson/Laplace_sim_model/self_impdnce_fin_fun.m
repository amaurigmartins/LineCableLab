function [Zm,Pm_0] = self_impdnce_fin_fun(x,h1,r1,he,omega,m0)

S12=2*(h1+he);
d12=r1;

P12_term1=x.*log((1+sqrt(1+(d12./x).^2))./(1+sqrt(1+(S12./x).^2)));
P12_term2=x.*log(S12/d12);
P12_term3=sqrt(x.^2+d12.^2);
P12_term4=sqrt(x.^2+S12.^2);
P12_term5=d12-S12;
Pm_0=P12_term1+P12_term2-P12_term3+P12_term4+P12_term5;

Zm=1i.*omega*m0./(2*pi).*Pm_0;