function [Pg_self]=P_papad_slf_2lay_under(h,r,e_g,m_g,sigma_g,t,f,con,kx)

if nargin==8
    kx=0; %sets default to zero
end

% Function for the 2-layered Self Earth Potential as per Papadopoulos' formula
% doi: 10.1049/iet-gtd.2010.0228

% Inputs
% h       : height of conductor [m]
% r       : external  radius of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% f       : frequency [rad/s]
% t       : layer thickness [m]
% Output
% Pg_self : Self earth potential coefficient of the underground conductor [Ohm/m]

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
e_0=e_g(1);
e_1=e_g(2); 
sigma_1=sigma_g(2);
e_2=e_g(3);
a_0=@(lambda) sqrt(lambda.^2+gamma_0.^2+k_x.^2);
a_1=@(lambda) sqrt(lambda.^2+gamma_1.^2+k_x.^2);
a_2=@(lambda) sqrt(lambda.^2+gamma_2.^2+k_x.^2);

% Smn = mu_n*a_m + mu_m*a_n
S10= @(lambda) m_0*a_1(lambda)+m_1*a_0(lambda);
S21= @(lambda) m_1*a_2(lambda)+m_2*a_1(lambda);

% Dmn = mu_m*a_n - mu_n*a_m
D10= @(lambda) m_1*a_0(lambda)-m_0*a_1(lambda);
D21= @(lambda) m_2*a_1(lambda)-m_1*a_2(lambda);

% Amn = @(lambda) a_n(lambda) * gamma_m^2 * mu_n + a_m(lambda) *gamma_n^2 * mu_m
A10 = @(lambda) a_0(lambda) * gamma_1^2 * m_0 + a_1(lambda) * gamma_0^2 * m_1;
A12 = @(lambda) a_2(lambda) * gamma_1^2 * m_2 + a_1(lambda) * gamma_2^2 * m_1;

% Deltamn = @(lambda) a_n(lambda) * gamma_m^2 * mu_n - a_m(lambda) * gamma_n^2 * mu_m
Delta10 = @(lambda) a_0(lambda) * gamma_1^2 * m_0 - a_1(lambda) * gamma_0^2 * m_1;
Delta12 = @(lambda) a_2(lambda) * gamma_1^2 * m_2 - a_1(lambda) * gamma_2^2 * m_1;


TOL=1e-3;

Pg_self=zeros(con,con);

for k=1:1:con
    h1=abs(h(1,k));
    h2=h1+TOL;
    dd=r(k);
    t=abs(t);

    % Self potential
    F=@(lambda) ( ...
            (S10(lambda)*S21(lambda)*exp(-a_1(lambda)*abs(h1-h2)) ...
            + S10(lambda)*D21(lambda)*exp(-a_1(lambda)*(2*t-h1-h2)) ...
            - D10(lambda)*S21(lambda)*exp(-a_1(lambda)*(h1+h2)) ...
            - D10(lambda)*D21(lambda)*exp(-a_1(lambda)*(2*t-abs(h1-h2)))) ...
            /...
            (a_1(lambda)*(S10(lambda)*S21(lambda)+D10(lambda)*D21(lambda)*exp(-2*a_1(lambda)*t))) ...
            );
    
    G1=@(lambda) ( ...
        (m_1*m_2*(gamma_1^2-gamma_2^2)*(S10(lambda)*A10(lambda)*exp(-a_1(lambda)*(2*t-h1-h2))-D10(lambda)*A10(lambda)*exp(-a_1(lambda)*(2*t+h1-h2))))...
        /...
        ((A10(lambda)*A12(lambda)-Delta10(lambda)*Delta12(lambda)*exp(-2*a_1(lambda)*t))*(S10(lambda)*S21(lambda)+D10(lambda)*D21(lambda)*exp(-2*a_1(lambda)*t)))...
        );

    G2=@(lambda) ( ...
        (m_1*m_2*(gamma_1^2-gamma_2^2)*(S10(lambda)*Delta10(lambda)*exp(-a_1(lambda)*(2*t+h2-h1))-D10(lambda)*Delta10(lambda)*exp(-a_1(lambda)*(2*t+h1+h2))))...
        /...
        ((A10(lambda)*A12(lambda)-Delta10(lambda)*Delta12(lambda)*exp(-2*a_1(lambda)*t))*(S10(lambda)*S21(lambda)+D10(lambda)*D21(lambda)*exp(-2*a_1(lambda)*t)))...
        );

    G3=@(lambda) ( ...
        (m_1*m_0*(gamma_1^2-gamma_0^2)*(S21(lambda)*Delta12(lambda)*exp(-a_1(lambda)*(2*t+h1-h2))+D21(lambda)*Delta12(lambda)*exp(-a_1(lambda)*(4*t-h1-h2))))...
        /...
        ((A10(lambda)*A12(lambda)-Delta10(lambda)*Delta12(lambda)*exp(-2*a_1(lambda)*t))*(S10(lambda)*S21(lambda)+D10(lambda)*D21(lambda)*exp(-2*a_1(lambda)*t)))...
        );

    G4=@(lambda) ( ...
        (m_1*m_0*(gamma_1^2-gamma_0^2)*(S21(lambda)*A12(lambda)*exp(-a_1(lambda)*(h1+h2))+D21(lambda)*A12(lambda)*exp(-a_1(lambda)*(2*t+h2-h1))))...
        /...
        ((A10(lambda)*A12(lambda)-Delta10(lambda)*Delta12(lambda)*exp(-2*a_1(lambda)*t))*(S10(lambda)*S21(lambda)+D10(lambda)*D21(lambda)*exp(-2*a_1(lambda)*t)))...
        );
    
    G=@(lambda) 2*a_1(lambda)*(G1(lambda)+G2(lambda)+G3(lambda)+G4(lambda));

    yfun=@(lambda) sum([0 (F(lambda)+G(lambda))*cos(dd*lambda)],'omitnan');
    Qs=integral(yfun,0,Inf,'ArrayValued',true);
    Pg_self(k,k)=(1i*w/(2*pi*(sigma_1+1i*w*e_1)))*Qs;
end
