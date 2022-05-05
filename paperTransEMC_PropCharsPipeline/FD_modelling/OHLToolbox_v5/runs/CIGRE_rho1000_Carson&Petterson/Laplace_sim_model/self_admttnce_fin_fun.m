function [P0] = self_admttnce_fin_fun(x,h1,r1)

D12=2*h1;
d12=r1;

P12_term1=x.*log((1+sqrt(1+(d12./x).^2))./(1+sqrt(1+(D12./x).^2)));
P12_term2=x.*log(D12/d12);
P12_term3=sqrt(x.^2+d12.^2);
P12_term4=sqrt(x.^2+D12.^2);
P12_term5=d12-D12;
P0=P12_term1+P12_term2-P12_term3+P12_term4+P12_term5;

