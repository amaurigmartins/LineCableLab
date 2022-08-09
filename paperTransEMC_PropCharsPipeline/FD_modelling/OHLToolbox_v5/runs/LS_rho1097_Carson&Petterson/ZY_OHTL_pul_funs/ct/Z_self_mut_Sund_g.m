%Imperfect ground by Sunde

function [Zg_Sund_mat]=Z_self_mut_Sund_g(h1,h2,h3,d12,d13,d23,e_g,m_g,sigma_g,omega)

% This function Calculates the Self and Mutual Earth Impedance by Sunde's
% Approxiamtion


% input: h1,h2,h3,e_g,_m_g,sigma_g,omega
% output: Z_self

m0=4*pi*1e-7;
gamma_g=sqrt(j*omega.*m_g.*(sigma_g+j*omega.*e_g));

% Self Impedance
ar_self=1+gamma_g.*h1;
par_self=gamma_g.*h1;

Zg11_self_Sund=j*omega.*m0/(2*pi).*log(ar_self./par_self);

ar_self=1+gamma_g.*h2;
par_self=gamma_g.*h2;

Zg22_self_Sund=j*omega.*m0/(2*pi).*log(ar_self./par_self);

ar_self=1+gamma_g.*h3;
par_self=gamma_g.*h3;

Zg33_self_Sund=j*omega.*m0/(2*pi).*log(ar_self./par_self);


% Mutual Impedance

ar_mut=(1+gamma_g.*((h1+h2)/2)).^2+(gamma_g.*(d12/2)).^2;
par_mut=(gamma_g.*((h1+h2)/2)).^2+(gamma_g.*(d12/2)).^2;

Zg12_mutual_Sund=j*omega.*m0/(4*pi).*log(ar_mut./par_mut);

ar_mut=(1+gamma_g.*((h2+h3)/2)).^2+(gamma_g.*(d23/2)).^2;
par_mut=(gamma_g.*((h2+h3)/2)).^2+(gamma_g.*(d23/2)).^2;

Zg23_mutual_Sund=j*omega.*m0/(4*pi).*log(ar_mut./par_mut);

ar_mut=(1+gamma_g.*((h1+h3)/2)).^2+(gamma_g.*(d13/2)).^2;
par_mut=(gamma_g.*((h1+h3)/2)).^2+(gamma_g.*(d13/2)).^2;

Zg13_mutual_Sund=j*omega.*m0/(4*pi).*log(ar_mut./par_mut);


Zg_Sund_mat=[Zg11_self_Sund Zg12_mutual_Sund Zg13_mutual_Sund; Zg12_mutual_Sund Zg22_self_Sund Zg23_mutual_Sund; Zg13_mutual_Sund Zg23_mutual_Sund Zg33_self_Sund];


