function [Zg_mutual]=Z_pol_mut_und(h,d,sigma1,f,con)

% Function for the Mutual Earth Impedance with underground-underground
% arrangement using Xue's formula

% Inputs
% h       : depth of conductor [m]
% d       : distances between conductors [m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]

% Output
% out: Mutual earth impedance between underground/underground conductors [Ohm/m]

% Constants
% sig0=0;
% eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;
ke=sqrt(1i*w*mu0*sigma1);

Zg_mutual=zeros(con,con);

if con == 1; return; end

% Mutual Impedance
for x=1:con
    for y=x+1:con
        if x~=y
            if h(1,x)<0 && h(1,y)<0
                h1=abs(h(1,x));
                h2=abs(h(1,y));
                dd=d(x,y);
                dnm=sqrt((h1-h2)^2+dd^2);
                Dnm=sqrt((h1+h2)^2+dd^2);
                K=(besselk(0,1i*ke*dnm,1)/exp(1i*ke*dnm))-(besselk(0,1i*ke*Dnm,1)/exp(1i*ke*Dnm));

                kernel=@(lambda) ((exp(-(h1+h2)*sqrt(lambda^2+ke^2)))/(lambda+sqrt(lambda^2+ke^2)))*cos(lambda*dd);
                % zfun=@(lambda) sum([0 zz(lambda)],'omitnan');
                Jm=integral(kernel,0,Inf,'ArrayValued',true);
                Zg_mutual(x,y)=(1i*w*mu0/(2*pi))*(K+2*Jm);
                Zg_mutual(y,x)=Zg_mutual(x,y);
            end
        end
    end
end
end

