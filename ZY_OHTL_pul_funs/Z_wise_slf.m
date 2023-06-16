% Z_wise_slf.m
function [Zg_self]=Z_wise_slf(h,e_g,m_g,sigma_g,omega,con)
% Function for the Self Earth Impedance by Wise's Integral
% Formula

% Inputs
% h1      : height of conductor [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Zg_self : Self earth impedance [Ohm/m]

% Constants
m0=4*pi*1e-7;
e0=8.854187817e-12;
k0=omega*(sqrt(m0*e0));
gama_g = sqrt(1i*omega*m0*(sigma_g + 1i*omega*e_g));

Zg_self=zeros(con,con);



for k=1:1:con
    if h(1,k) >0
    % Self Impedance
    H=2*h(1,k);
    z_lambda =@(lambda) sum([0 ((2*m_g*(exp(-H*lambda)))/((sqrt(lambda^2 + gama_g^2 + k0^2))*m0 + lambda*m_g))],'omitnan');
    Jm = integral(z_lambda,0,Inf,'ArrayValued',true);

    Zg_self(k,k)=1i*omega*m0/(2*pi)*(+Jm); % let us add the perfect earth term later, to keep consistency

    end

end