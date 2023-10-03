function [Pg_self]=P_ohl_slf_2lay(h,r,e_g,m_g,sigma_g,t,f,con,kx)

if nargin==8
    kx=0; %sets default to zero
end

% Function for the 2-layered Self Earth Potential as per Sunde's formula

% Inputs
% h       : height of conductor [m]
% r       : external  radius of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% f       : frequency [rad/s]
% t       : layer thickness [m]
% Output
% Pg_self : Self earth potential coefficient of the overhead conductor [Ohm/m]

% Constants
sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;

e_g = [eps0 e_g];
m_g = [mu0 m_g];
sigma_g = [sig0 sigma_g];

if kx==0
    k_x=@(omega) omega.*0;  
elseif strcmp(kx,'k0')
    k_x=@(omega) omega.*sqrt(m_g(1).*(e_g(1)-1i.*(sigma_g(1)./omega)));
elseif strcmp(kx,'k1')
    k_x=@(omega) omega.*sqrt(m_g(2).*(e_g(2)-1i.*(sigma_g(2)./omega)));
elseif strcmp(kx,'k2')
    k_x=@(omega) omega.*sqrt(m_g(3).*(e_g(3)-1i.*(sigma_g(3)./omega)));
end

gamma_1=@(omega) sqrt(1i.*omega.*m_g(1).*(sigma_g(1)+1i.*omega.*e_g(1)));
gamma_2=@(omega) sqrt(1i.*omega.*m_g(2).*(sigma_g(2)+1i.*omega.*e_g(2)));
gamma_3=@(omega) sqrt(1i.*omega.*m_g(3).*(sigma_g(3)+1i.*omega.*e_g(3)));
a_1=@(lambda,omega) sqrt(lambda.^2+gamma_1(omega).^2+k_x(omega).^2);
a_2=@(lambda,omega) sqrt(lambda.^2+gamma_2(omega).^2+k_x(omega).^2);
a_3=@(lambda,omega) sqrt(lambda.^2+gamma_3(omega).^2+k_x(omega).^2);

TOL=1e-3;

qq1= @(a1,a2,a3,gamma1,gamma2,gamma3,hi,mu1,mu2,mu3,omega,t2)(1.0./gamma1.^2.*mu1.^2.*mu2.*omega.*exp(a2.*t2).*(a2.*a3.*gamma2.^4.*mu3.^2.*exp(a1.*hi).*1i+a3.^2.*gamma2.^4.*mu2.*mu3.*exp(a1.*hi).*1i-a2.*a3.*gamma2.^4.*mu3.^2.*exp(a1.*hi+a2.*t2.*4.0).*1i-a3.^2.*gamma2.^4.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*2.0).*2.0i+a3.^2.*gamma2.^4.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*4.0).*1i-a2.*a3.*gamma1.^2.*gamma2.^2.*mu3.^2.*exp(a1.*hi).*1i-a2.*a3.*gamma1.^2.*gamma3.^2.*mu2.^2.*exp(a1.*hi).*1i+a2.*a3.*gamma2.^2.*gamma3.^2.*mu2.^2.*exp(a1.*hi).*1i-a2.^2.*gamma1.^2.*gamma3.^2.*mu2.*mu3.*exp(a1.*hi).*1i-a3.^2.*gamma1.^2.*gamma2.^2.*mu2.*mu3.*exp(a1.*hi).*1i+a2.^2.*gamma2.^2.*gamma3.^2.*mu2.*mu3.*exp(a1.*hi).*1i+a2.*a3.*gamma1.^2.*gamma2.^2.*mu3.^2.*exp(a1.*hi+a2.*t2.*4.0).*1i+a2.*a3.*gamma1.^2.*gamma3.^2.*mu2.^2.*exp(a1.*hi+a2.*t2.*4.0).*1i-a2.*a3.*gamma2.^2.*gamma3.^2.*mu2.^2.*exp(a1.*hi+a2.*t2.*4.0).*1i-a2.^2.*gamma1.^2.*gamma2.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*2.0).*4.0i+a2.^2.*gamma1.^2.*gamma3.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*2.0).*2.0i+a3.^2.*gamma1.^2.*gamma2.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*2.0).*2.0i+a2.^2.*gamma2.^2.*gamma3.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*2.0).*2.0i-a2.^2.*gamma1.^2.*gamma3.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*4.0).*1i-a3.^2.*gamma1.^2.*gamma2.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*4.0).*1i+a2.^2.*gamma2.^2.*gamma3.^2.*mu2.*mu3.*exp(a1.*hi+a2.*t2.*4.0).*1i).*(-1.0./2.0))./(pi.*(a1.*a3.*mu2.^2.*exp(a1.*hi.*2.0+a2.*t2)-a1.*a3.*mu2.^2.*exp(a1.*hi.*2.0+a2.*t2.*3.0)+a2.^2.*mu1.*mu3.*exp(a1.*hi.*2.0+a2.*t2)-a2.^2.*mu1.*mu3.*exp(a1.*hi.*2.0+a2.*t2.*3.0)+a1.*a2.*mu2.*mu3.*exp(a1.*hi.*2.0+a2.*t2)+a2.*a3.*mu1.*mu2.*exp(a1.*hi.*2.0+a2.*t2)+a1.*a2.*mu2.*mu3.*exp(a1.*hi.*2.0+a2.*t2.*3.0)+a2.*a3.*mu1.*mu2.*exp(a1.*hi.*2.0+a2.*t2.*3.0)).*(a2.^2.*gamma1.^2.*gamma3.^2.*mu2.^2-a2.^2.*gamma1.^2.*gamma3.^2.*mu2.^2.*exp(a2.*t2.*2.0)+a1.*a3.*gamma2.^4.*mu1.*mu3-a1.*a3.*gamma2.^4.*mu1.*mu3.*exp(a2.*t2.*2.0)+a1.*a2.*gamma2.^2.*gamma3.^2.*mu1.*mu2+a2.*a3.*gamma1.^2.*gamma2.^2.*mu2.*mu3+a1.*a2.*gamma2.^2.*gamma3.^2.*mu1.*mu2.*exp(a2.*t2.*2.0)+a2.*a3.*gamma1.^2.*gamma2.^2.*mu2.*mu3.*exp(a2.*t2.*2.0)));

gg1=@(a1,a2,a3,gamma1,hi,mu1,mu2,mu3,omega,t2)(1.0./gamma1.^2.*mu1.*omega.*exp(-a1.*hi).*(a1.*a3.*mu2.^2-a2.^2.*mu1.*mu3-a1.*a3.*mu2.^2.*exp(a2.*t2.*2.0)+a2.^2.*mu1.*mu3.*exp(a2.*t2.*2.0)+a1.*a2.*mu2.*mu3-a2.*a3.*mu1.*mu2+a1.*a2.*mu2.*mu3.*exp(a2.*t2.*2.0)-a2.*a3.*mu1.*mu2.*exp(a2.*t2.*2.0)).*2.5e-1i)./(a1.*pi.*(a1.*a3.*mu2.^2+a2.^2.*mu1.*mu3-a1.*a3.*mu2.^2.*exp(a2.*t2.*2.0)-a2.^2.*mu1.*mu3.*exp(a2.*t2.*2.0)+a1.*a2.*mu2.*mu3+a2.*a3.*mu1.*mu2+a1.*a2.*mu2.*mu3.*exp(a2.*t2.*2.0)+a2.*a3.*mu1.*mu2.*exp(a2.*t2.*2.0)));

yy=@(a1,g1,gamma1,hi,hj,lambda,mu1,omega,q1,y)cos(lambda.*y).*(g1.*exp(-a1.*hj).*2.0+a1.*q1.*exp(-a1.*hj).*2.0+(1.0./gamma1.^2.*mu1.*omega.*exp(-a1.*abs(hi-hj+TOL)).*5.0e-1i)./(a1.*pi));


Pg_self=zeros(con,con);

for k=1:1:con

    if h(1,k) > 0

        % Self potential
        yfun=@(lambda,omega) sum([0 ...
            yy(a_1(lambda,omega),...
            gg1(a_1(lambda,omega),a_2(lambda,omega),a_3(lambda,omega),gamma_1(omega),h(1,k),m_g(1),m_g(2),m_g(3),omega,t),...
            gamma_1(omega),...
            h(1,k),...
            h(1,k),...
            lambda,...
            m_g(1),...
            omega,...
            qq1(a_1(lambda,omega),a_2(lambda,omega),a_3(lambda,omega),gamma_1(omega),gamma_2(omega),gamma_3(omega),h(1,k),m_g(1),m_g(2),m_g(3),omega,t),...
            r(k))...
            ],'omitnan');
        yfun=@(lambda) yfun(lambda,w);
        Qs=integral(yfun,0,Inf,'ArrayValued',true);
        Pg_self(k,k)=1i*w*Qs;
    end

end