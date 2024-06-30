function [Zg_mutual]=Z_papad_mut_2lay_under(h,d,e_g,m_g,sigma_g,t,f,con,kx)

if nargin==8
    kx=0; %sets default to zero
end

% Function for the 2-layered Mutual Earth Impedance as per Papadopoulos' formula
% doi: 10.1049/iet-gtd.2010.0228

% Inputs
% h       : height of conductor [m]
% D       : distances between conductors [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% f       : frequency [rad/s]
% t       : layer thickness [m]
% Output
% Zg_mutual : Mutual earth impedance of the underground conductor [Ohm/m]

% Constants
sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;

e_g = [eps0; e_g];
m_g = [mu0; m_g];
sigma_g = [sig0; sigma_g];

if kx==0
    k_x= w.*0;
else
    k_x= w.*sqrt(m_g(kx).*(e_g(kx)-1i.*(sigma_g(kx)./w)));
end

gamma_0= sqrt(1i.*w.*m_g(1).*(sigma_g(1)+1i.*w.*e_g(1)));
gamma_1= sqrt(1i.*w.*m_g(2).*(sigma_g(2)+1i.*w.*e_g(2)));
gamma_2= sqrt(1i.*w.*m_g(3).*(sigma_g(3)+1i.*w.*e_g(3)));
m_0=m_g(1);
m_1=m_g(2);
m_2=m_g(3);
a_0=@(lambda) sqrt(lambda.^2+gamma_0.^2+k_x.^2);
a_1=@(lambda) sqrt(lambda.^2+gamma_1.^2+k_x.^2);
a_2=@(lambda) sqrt(lambda.^2+gamma_2.^2+k_x.^2);

% Smn = mu_n*a_m + mu_m*a_n
% Dmn = mu_m*a_n - mu_n*a_m
% Amn = @(lambda) a_n(lambda) * gamma_m^2 * mu_n + a_m(lambda) *gamma_n^2 * mu_m
% Deltamn = @(lambda) a_n(lambda) * gamma_m^2 * mu_n - a_m(lambda) * gamma_n^2 * mu_m
S10= @(lambda) m_0*a_1(lambda)+m_1*a_0(lambda);
S21= @(lambda) m_1*a_2(lambda)+m_2*a_1(lambda);
D10= @(lambda) m_1*a_0(lambda)-m_0*a_1(lambda);
D21= @(lambda) m_2*a_1(lambda)-m_1*a_2(lambda);

TOL=1e-3;

Zg_mutual=zeros(con,con);

% Mutual Impedance
for x=1:1:con
    for y=1:1:con
        if x~=y
            h1=abs(h(1,x));
            h2=abs(h(1,y));
            dd=d(x,y);
            t=abs(t);
            if abs(h2-h1)<TOL; h2=h2+TOL;end
            F=@(lambda) ( ...
                (S10(lambda)*S21(lambda)*exp(-a_1(lambda)*abs(h1-h2)) ...
                + S10(lambda)*D21(lambda)*exp(-a_1(lambda)*(2*t-h1-h2)) ...
                - D10(lambda)*S21(lambda)*exp(-a_1(lambda)*(h1+h2)) ...
                - D10(lambda)*D21(lambda)*exp(-a_1(lambda)*(2*t-abs(h1-h2)))) ...
                /...
                (a_1(lambda)*(S10(lambda)*S21(lambda)+D10(lambda)*D21(lambda)*exp(-2*a_1(lambda)*t))) ...
                );
            zfun=@(lambda) sum([0 F(lambda)*cos(dd*lambda)],'omitnan');
            Jm=integral(zfun,0,Inf,'ArrayValued',true);

            Zg_mutual(x,y)=((1i*w*m_1)/(2*pi))*Jm;
        end
    end
end


end