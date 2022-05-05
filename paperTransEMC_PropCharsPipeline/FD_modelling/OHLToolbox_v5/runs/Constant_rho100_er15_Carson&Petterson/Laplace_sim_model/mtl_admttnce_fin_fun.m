function [P0] = mtl_admttnce_fin_fun(x,h1,h2,y)

D12=sqrt((h1+h2).^2+y.^2);
d12=sqrt((h1-h2).^2+y.^2);

P12_term1=x.*log((1+sqrt(1+(d12./x).^2))./(1+sqrt(1+(D12./x).^2)));
P12_term2=x.*log(D12/d12);
P12_term3=sqrt(x.^2+d12.^2);
P12_term4=sqrt(x.^2+D12.^2);
P12_term5=d12-D12;
P0=P12_term1+P12_term2-P12_term3+P12_term4+P12_term5;

% C=2*pi*e0.*(x.^2./P0);
% Y=1i.*omega*C;

