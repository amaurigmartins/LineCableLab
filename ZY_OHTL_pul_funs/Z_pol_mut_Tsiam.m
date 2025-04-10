function [Zg_mutual]=Z_pol_mut_Tsiam(h,d,eps1,mu1,sigma1,f,con)

% Function for the Mutual Earth Impedance with overhead-underground
% arrangement using Pollaczek formula (Tsiamitros eq. 45)

% Inputs
% h       : depth of conductor [m]
% d       : distances between conductors [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% f       : frequency [Hz]

% Output
% out: Mutual earth impedance between overhead/underground conductors [Ohm/m]

% Constants
sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;
k_x=@(omega) 0;

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

            if (h1 < 0 && h2 > 0) %Z10
                zz=@(a0,a1,hi,hj,lambda,mu0,mu1,omega,y) ...
                    (mu0.*mu1.*omega.*exp(a1.*hi-a0.*hj).*cos(lambda.*y) ...
                    .*1i)./(pi.*(a0.*mu1+a1.*mu0));
            elseif (h1 > 0 && h2 < 0) %Z01
                zz=@(a0,a1,hi,hj,lambda,mu0,mu1,omega,y) ...
                    (mu0.*mu1.*omega.*exp(-a0.*hi+a1.*hj).*cos(lambda.*y) ...
                    .*1i)./(pi.*(a0.*mu1+a1.*mu0));
            else
                continue
            end

            zfun=@(lambda,omega) sum([0 zz(a_0(lambda,omega),a_1(lambda,omega),h1,h2,lambda,mu0,mu1,omega,d(x,y))],'omitnan');
            zfun=@(lambda) zfun(lambda,w);
            Jm=integral(zfun,0,Inf,'ArrayValued',true);

            Zg_mutual(x,y)=Jm;
            Zg_mutual(y,x)=Zg_mutual(x,y);

        end
    end
end
end

