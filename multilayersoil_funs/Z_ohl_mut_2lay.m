function [Zg_mutual]=Z_ohl_mut_2lay(h,d,e_g,m_g,sigma_g,t,f,con,kx)

if nargin==8
    kx=0; %sets default to zero
end

% Function for the 2-layered Mutual Earth Impedance as per Sunde's formula

% Inputs
% h       : height of conductor [m]
% D       : distances between conductors [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% f       : frequency [rad/s]
% t       : layer thickness [m]
% Output
% Zg_mutual : Mutual earth impedance of the overhead conductor [Ohm/m]

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

gg1=@(a1,a2,a3,gamma1,hi,mu1,mu2,mu3,omega,t2)...
    (1.0./gamma1.^2.*mu1.*omega.*exp(-a1.*hi).*...
    (a1.*a3.*mu2.^2-a2.^2.*mu1.*mu3-a1.*a3.*mu2.^2.*exp(a2.*t2.*2.0)+...
    a2.^2.*mu1.*mu3.*exp(a2.*t2.*2.0)+a1.*a2.*mu2.*mu3-a2.*a3.*mu1.*mu2+...
    a1.*a2.*mu2.*mu3.*exp(a2.*t2.*2.0)-a2.*a3.*mu1.*mu2.*exp(a2.*t2.*2.0)).*2.5e-1i)./...
    (a1.*pi.*(a1.*a3.*mu2.^2+a2.^2.*mu1.*mu3-a1.*a3.*mu2.^2.*exp(a2.*t2.*2.0)-...
    a2.^2.*mu1.*mu3.*exp(a2.*t2.*2.0)+a1.*a2.*mu2.*mu3+a2.*a3.*mu1.*mu2+...
    a1.*a2.*mu2.*mu3.*exp(a2.*t2.*2.0)+a2.*a3.*mu1.*mu2.*exp(a2.*t2.*2.0)));

zz=@(a1,g1,gamma1,hi,hj,lambda,mu1,omega,y)...
    gamma1.^2.*cos(lambda.*y).*(g1.*exp(-a1.*hj)+...
    (1.0./gamma1.^2.*mu1.*omega.*exp(-a1.*abs(hi-hj+TOL)).*2.5e-1i)./(a1.*pi)).*2.0;


Zg_mutual=zeros(con,con);

% Mutual Impedance
for x=1:1:con
    for y=1:1:con
        if x~=y
            h1=h(1,x);
            h2=h(1,y);

            if (h1 > 0 && h2 > 0)

                zfun=@(lambda,omega) sum([0 ...
                    zz(a_1(lambda,omega),...
                    gg1(...
                        a_1(lambda,omega),...
                        a_2(lambda,omega),...
                        a_3(lambda,omega),...
                        gamma_1(omega),...
                        h1,...
                        m_g(1),...
                        m_g(2),...
                        m_g(3),...
                        omega,...
                        t),...
                    gamma_1(omega),...
                    h1,...
                    h2,...
                    lambda,...
                    m_g(1),...
                    omega,...
                    d(x,y))...
                    ],'omitnan');
                zfun=@(lambda) zfun(lambda,w);
                Jm=integral(zfun,0,Inf,'ArrayValued',true);

                Zg_mutual(x,y)=Jm;

            end
        end
    end
end


end