% Z_wise_mut.m
function [Zg_mutual]=Z_wise_mut(h,d,e_g,m_g,sigma_g,omega,con)
% Function for the Mutual Earth Impedance Impedance by Wise's Integral
% Formula

% Inputs
% h1        : height of conductor [m]
% h2        : height of conductor [m]
% d         : separation distance [m]
% e_g       : permittivity of earth [F/m]
% m_g       : permeability of earth [H/m]
% sigma_g   : conductivity of earth [S/m]
% omega     : angular frequency [rad/s]
% Output
% Zg_mutual : Mutual earth impedance [Ohm/m]


m0=4*pi*1e-7;
e0=8.854187817e-12;
k0=omega*(sqrt(m0*e0));
gama_g = sqrt(1i*omega*m0*(sigma_g + 1i*omega*e_g));

Zg_mutual=zeros(con,con);

% Mutual Impedance
for x=1:con
    for y=x+1:con
        if x~=y
            h1=h(1,x);
            h2=h(1,y);
            H=h1+h2;
            
            z_lambda =@(lambda) sum([0 ((2*m_g*(exp(-H*lambda)))/((sqrt(lambda^2 + gama_g^2 + k0^2))*m0 + lambda*m_g))*cos(lambda*d(x,y))], 'omitnan');
            Jm = integral(z_lambda,0,Inf,'ArrayValued',true);

            Zg_mutual(x,y)=1i*omega*m0/(2*pi)*(+Jm);
            Zg_mutual(y,x)=Zg_mutual(x,y);
        end
    end
end

