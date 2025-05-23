function [Zg_mutual]=Z_papad_mut(h,d,eps1,mu1,sigma1,f,con,kx)

if nargin==7
    kx=0; %sets default to zero
end

% Function for the Mutual Earth Impedance as per Papadopoulos formula

% Inputs
% h       : depth of conductor [m]
% d       : distances between conductors [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]

% Output
% out: Mutual earth impedance between underground conductors [Ohm/m]

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

Zg_mutual=zeros(con,con);

if con == 1; return; end

% Mutual Impedance
for x=1:con
    for y=x+1:con
        if x~=y
            h1=h(1,x);
            h2=h(1,y);

            if (h1 < 0 && h2 <0)

                zz=@(a0,a1,hi,hj,lambda,mu0,mu1,omega,y) ...
                    (mu1.*omega.*exp(-a1.*abs(hi-hj+1e-3)).*cos(lambda.*y).*5.0e-1i)./(a1.*pi) ...
                    -(mu1.*omega.*exp(a1.*(hi+hj)).*cos(lambda.*y).*(a0.*mu1+a1.*mu0.*sign(hi)).*5.0e-1i) ...
                    ./(a1.*pi.*(a0.*mu1+a1.*mu0));


                zfun=@(lambda,omega) sum([0 zz(a_0(lambda,omega),a_1(lambda,omega),h1,h2,lambda,mu0,mu1,omega,d(x,y))],'omitnan');
                zfun=@(lambda) zfun(lambda,w);
                Jm=integral(zfun,0,Inf,'ArrayValued',true);

                Zg_mutual(x,y)=Jm;
                Zg_mutual(y,x)=Zg_mutual(x,y);

            end
        end
    end
end



end
