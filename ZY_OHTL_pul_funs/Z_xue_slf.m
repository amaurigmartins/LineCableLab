function [Zg_self]=Z_xue_slf(h,r,eps1,sigma1,f,con)

% Function for the Self Earth Impedance as per Xue's formula

% Inputs
% h       : depth of conductor [m]
% r       : external  radius of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]

% Output
% out: Self earth impedance of the underground conductor [Ohm/m]

% Constants
% sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;
ke=sqrt(-1i*w*mu0*(sigma1+1i*w*eps1));
ka=sqrt(w^2*mu0*eps0);
% kv=0;
u0=@(lambda) sqrt(lambda^2-ka^2);
u1=@(lambda) sqrt(lambda^2-ke^2);

Zg_self=zeros(con,con);

for k=1:1:con
    if h(1,k) < 0
        h1=h(1,k);
        h2=h(1,k);
        dd=r(k);
        dnm=sqrt((h1-h2)^2+dd^2);
        Dnm=sqrt((h1+h2)^2+dd^2);
        K=besselk(0,1i*ke*dnm)-besselk(0,1i*ke*Dnm);
        T4=@(lambda) (exp((h1+h2)*u1(lambda))/(u0(lambda)+u1(lambda)))*(lambda^2/u1(lambda)^2);
        T6=@(lambda) (exp((h1+h2)*u1(lambda))/(u0(lambda)+u1(lambda)))*(1/u1(lambda)^2);
        zz=@(lambda) (2*T4(lambda)-2*ke^2*T6(lambda))*cos(lambda*dd);
        zfun=@(lambda) sum([0 zz(lambda)],'omitnan');
        Js=integral(zfun,0,Inf,'ArrayValued',true);
        Zg_self(k,k)=(1i*w*mu0/(2*pi))*(K+Js);
    end

end