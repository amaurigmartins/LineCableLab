function [Pg_self]=P_papad_slf(h,r,eps1,mu1,sigma1,f,con,kx)

if nargin==7
    kx=0; %sets default to zero
end

% Function for the Self Earth Potential Coefficient as per Papadopoulos formula

% Inputs
% h       : depth of conductor [m]
% r       : external  radius of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]

% Output
% out: Self earth potential of the underground conductor [Ohm/m]

% Constants
sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;

if strcmp(kx,'k1')
    k_x=@(omega) omega.*sqrt(mu1.*(eps1-1i.*(sigma1./omega)));
elseif strcmp(kx,'k0')
    k_x=@(omega) omega.*sqrt(mu0.*eps0);
else
    k_x=@(omega) omega.*0;
end

gamma_0=@(omega) sqrt(1i.*omega.*mu0.*(sig0+1i.*omega.*eps0));
gamma_1=@(omega) sqrt(1i.*omega.*mu1.*(sigma1+1i.*omega.*eps1));
a_0=@(lambda,omega) sqrt(lambda.^2+gamma_0(omega).^2+k_x(omega).^2);
a_1=@(lambda,omega) sqrt(lambda.^2+gamma_1(omega).^2+k_x(omega).^2);

Pg_self=zeros(con,con);

for k=1:1:con

    if h(1,k) < 0
        % Self potential coefficient
        yy=@(a0,a1,gamma0,gamma1,hi,hj,lambda,mu0,mu1,omega,y) ...
            (1.0./gamma1.^2.*mu1.*omega.*exp(-a1.*abs(hi-hj+1e-3)).*cos(lambda.*y) ...
            .*5.0e-1i)./(a1.*pi)-(1.0./gamma1.^2.*mu1.*omega.*exp(a1.*(hi+hj)) ...
            .*cos(lambda.*y).*(a0.*mu1+a1.*mu0.*sign(hi)).*5.0e-1i) ...
            ./(a1.*pi.*(a0.*mu1+a1.*mu0))+(a1.*1.0./gamma1.^2.*mu0.*mu1.^2.*omega.*exp(a1.*(hi+hj)) ...
            .*cos(lambda.*y).*(sign(hi)-1.0).*(gamma0.^2-gamma1.^2).*5.0e-1i) ...
            ./(pi.*(a0.*gamma1.^2.*mu0+a1.*gamma0.^2.*mu1).*(a0.*mu1+a1.*mu0));

        yfun=@(lambda,omega) sum([0 yy(a_0(lambda,omega),a_1(lambda,omega),gamma_0(omega),gamma_1(omega),h(1,k),h(1,k),lambda,mu0,mu1,omega,r(k))],'omitnan');
        yfun=@(lambda) yfun(lambda,w);

        Qs=integral(yfun,0,Inf,'ArrayValued',true);
        Pg_self(k,k)=1i*w*Qs;
    end

end