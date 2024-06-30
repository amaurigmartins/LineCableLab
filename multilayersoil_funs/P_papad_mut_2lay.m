function [Pg_mutual]=P_papad_mut_2lay(h,d,e_g,m_g,sigma_g,t,f,con,kx)

if nargin==8
    kx=0; %sets default to zero
end

% Function for the 2-layered Mutual Earth Potential as per Papadopoulos' formula

% Inputs
% h       : height of conductor [m]
% D       : distances between conductors [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% f       : frequency [rad/s]
% t       : layer thickness [m]
% Output
% Pg_mutual : Mutual earth potential of the overhead conductor [Ohm/m]

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
e_2=e_g(3);

a_0=@(lambda) sqrt(lambda.^2+gamma_0.^2+k_x.^2);
a_1=@(lambda) sqrt(lambda.^2+gamma_1.^2+k_x.^2);
a_2=@(lambda) sqrt(lambda.^2+gamma_2.^2+k_x.^2);

s01=@(lambda) a_0(lambda)*m_1 + a_1(lambda)*m_0;
d01=@(lambda) a_0(lambda)*m_1 - a_1(lambda)*m_0;
s12=@(lambda) a_1(lambda)*m_2 + a_2(lambda)*m_1;
d12=@(lambda) a_1(lambda)*m_2 - a_2(lambda)*m_1;

S01=@(lambda) m_0*gamma_1^2*a_0(lambda) + m_1*gamma_0^2*a_1(lambda);
D01=@(lambda) m_0*gamma_1^2*a_0(lambda) - m_1*gamma_0^2*a_1(lambda);
S12=@(lambda) m_1*gamma_2^2*a_1(lambda) + m_2*gamma_1^2*a_2(lambda);
D12=@(lambda) m_1*gamma_2^2*a_1(lambda) - m_2*gamma_1^2*a_2(lambda);

Delta=@(lambda)  s01(lambda)*s12(lambda)+d01(lambda)*d12(lambda)*exp(-2*a_1(lambda)*abs(t));
Delta2=@(lambda) S01(lambda)*S12(lambda)+D01(lambda)*D12(lambda)*exp(-2*a_1(lambda)*abs(t));


TOL=1e-3;

Pg_mutual=zeros(con,con);

for x=1:1:con
    for y=1:1:con
        if x~=y
            h1=h(1,x);
            h2=h(1,y);
            dd=d(x,y);
            d1=sqrt(dd.^2+(h1-h2).^2);
            d2=sqrt(dd.^2+(h1+h2).^2); %log(d2./d1)
            

            if (h1 > 0 && h2 > 0)
                Fstrat=@(lambda) m_0*((s12(lambda)+d12(lambda)*exp(-2*a_1(lambda)*abs(t))) ...
                    /...
                    (s01(lambda)*s12(lambda)+d01(lambda)*d12(lambda)*exp(-2*a_1(lambda)*abs(t))));
                
                Gstrat=@(lambda) lambda*( ...
                    m_0*m_1*(gamma_0^2-gamma_1^2) * (s12(lambda)+d12(lambda)*exp(-2*a_1(lambda)*abs(t))) * (S12(lambda)+D12(lambda)*exp(-2*a_1(lambda)*abs(t))) ...
                    -4*m_0*m_1^2*m_2*a_1(lambda)^2*gamma_0^2*(gamma_2^2-gamma_1^2)*exp(-2*a_1(lambda)*abs(t)) ...
                    ) ./ (Delta(lambda)*Delta2(lambda));
                    
                
                yfun=@(lambda) sum([0 (Fstrat(lambda)+Gstrat(lambda))*exp(-lambda*(h1+h2))*cos(dd*lambda)],'omitnan');
                Qm=integral(yfun,0,Inf,'ArrayValued',true);

                % Pg_mutual(x,y)=1i*w*Qm;
                Pg_mutual(x,y)=((1/(2*pi*e_0))*log(d2./d1)+(1/(pi*e_0))*Qm);

            end
        end
    end
end


end