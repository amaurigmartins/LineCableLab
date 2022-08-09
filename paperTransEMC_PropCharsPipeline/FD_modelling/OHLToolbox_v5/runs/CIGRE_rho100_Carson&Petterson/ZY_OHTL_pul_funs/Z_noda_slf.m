function [Zg_self]=Z_noda_slf(h,e_g,sigma_g,omega,r_ext,con)
% Function for the Self Earth Impedance by Noda's Approximate
% Formula

% Inputs
% h1      : height of conductor [m]
% r_ext   : conductor radius    [m]
% e_g     : permittivity of earth [F/m] (originally not used)
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Zg_self : Self earth impedance [Ohm/m]

m0=4*pi*1e-7;
D1=r_ext;

Zg_self=zeros(con,con);

for k=1:1:con
    H=2*h(1,k);
    x=0;
    p=1/sqrt(1i*omega*m0.*sigma_g);
    theta=rad2deg(atan(x/H));
    if theta <= 50.45
        A=0.07360;
        a=0.1500;
    else
        A=0.00247*theta-0.05127;
        a=0.004726*theta-0.08852;
    end    

    % Self Impedance
    beta=(1-A*a)/(1-A);
    term1=A.*log(sqrt((H+2.*a.*p).^2+x^2)./D1(con));
    term2=(1-A).*log(sqrt((H+2.*beta.*p).^2+x^2)./D1(con));

    Zg_self(k,k)=1i*omega.*m0/(2*pi).*(term1+term2);
end