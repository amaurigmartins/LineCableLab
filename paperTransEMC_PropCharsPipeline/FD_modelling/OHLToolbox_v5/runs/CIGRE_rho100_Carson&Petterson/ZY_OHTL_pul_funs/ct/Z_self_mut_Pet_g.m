%Imperfect ground by Pettersson

function [Zg_Pet_mat]=Z_self_mut_Pet_g(h1,h2,h3,d12,d13,d23,e_g,m_g,sigma_g,omega)

% This function Calculates the Self and Mutual Earth Impedance by Pettersson's
% Approxiamtion

% input: h1,e_g,_m_g,sigma_g,omega
% output: Z_self
m0=4*pi*1e-7;
e0=8.854187817e-12;
b=sqrt(j*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));


% Self Impedance
ar_self=1+b*h1;
par_self=b*h1;

Zg11_self_Pet=j*omega.*m0./(2*pi).*log(ar_self./par_self);

ar_self=1+b*h2;
par_self=b*h2;

Zg22_self_Pet=j*omega.*m0./(2*pi).*log(ar_self./par_self);

ar_self=1+b*h3;
par_self=b*h3;

Zg33_self_Pet=j*omega.*m0./(2*pi).*log(ar_self./par_self);


% Mutual Impedance
dp=sqrt((h1+h2+2./b).^2+d12.^2);
d2=sqrt(d12.^2+(h1+h2).^2);
term=log(dp./d2);


Zg12_mutual_Pet=j*omega.*m0./(2*pi).*term;


dp=sqrt((h1+h3+2./b).^2+d13.^2);
d2=sqrt(d13.^2+(h1+h3).^2);
term=log(dp./d2);


Zg13_mutual_Pet=j*omega.*m0./(2*pi).*term;


dp=sqrt((h2+h3+2./b).^2+d23.^2);
d2=sqrt(d23.^2+(h2+h3).^2);
term=log(dp./d2);


Zg23_mutual_Pet=j*omega.*m0./(2*pi).*term;


Zg_Pet_mat=[Zg11_self_Pet Zg12_mutual_Pet Zg13_mutual_Pet; Zg12_mutual_Pet Zg22_self_Pet Zg23_mutual_Pet; Zg13_mutual_Pet Zg23_mutual_Pet Zg33_self_Pet];


