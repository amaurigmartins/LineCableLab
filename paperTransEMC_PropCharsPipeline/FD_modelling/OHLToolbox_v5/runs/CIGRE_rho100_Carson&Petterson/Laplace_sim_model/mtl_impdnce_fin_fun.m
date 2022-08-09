function [Zm,Pm_0] = mtl_impdnce_fin_fun(x,h1,h2,y,he,omega,m0)

S12=sqrt((h1+h2+2*he).^2+y.^2);
d12=sqrt((h1-h2).^2+y.^2);

P12_term1=x.*log((1+sqrt(1+(d12./x).^2))./(1+sqrt(1+(S12./x).^2)));
P12_term2=x.*log(S12/d12);
P12_term3=sqrt(x.^2+d12.^2);
P12_term4=sqrt(x.^2+S12.^2);
P12_term5=d12-S12;
Pm_0=P12_term1+P12_term2-P12_term3+P12_term4+P12_term5;

Zm=1i.*omega*m0./(2*pi).*Pm_0;