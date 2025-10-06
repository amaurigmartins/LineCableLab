function [Zg_self]=Z_pol_slf_und(h,r,sigma1,f,con)

% Function for the Self Earth Impedance as per Pollaczek's formula

% Inputs
% h       : depth of conductor [m]
% r       : external  radius of conductor [m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]

% Output
% out: Self earth impedance of the underground conductor [Ohm/m]

% Constants
% sig0=0;
% eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;
ke=sqrt(1i*w*mu0*sigma1);


Zg_self=zeros(con,con);

for k=1:1:con
    if h(1,k) < 0
        h1=abs(h(1,k));
        h2=abs(h(1,k));
        dd=r(k);
        dnm=sqrt((h1-h2)^2+dd^2);
        Dnm=sqrt((h1+h2)^2+dd^2);
        K=(besselk(0,1i*ke*dnm,1)/exp(1i*ke*dnm))-(besselk(0,1i*ke*Dnm,1)/exp(1i*ke*Dnm));
        
        kernel=@(lambda) ((exp(-(h1+h2)*sqrt(lambda^2+ke^2)))/(lambda+sqrt(lambda^2+ke^2)))*cos(lambda*dd);
        % zfun=@(lambda) sum([0 zz(lambda)],'omitnan');
        Js=integral(kernel,0,Inf,'ArrayValued',true);
        Zg_self(k,k)=(1i*w*mu0/(2*pi))*(K+2*Js);
        
    end

end