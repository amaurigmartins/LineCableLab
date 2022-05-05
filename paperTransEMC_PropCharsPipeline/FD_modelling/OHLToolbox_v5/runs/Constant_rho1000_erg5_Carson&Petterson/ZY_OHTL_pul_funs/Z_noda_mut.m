function [Zg_self]=Z_noda_mut(h,e_g,sigma_g,omega,d,con)
% Function for the Self Earth Impedance by Noda's Approximate
% Formula

% Inputs
% h1      : height of conductor [m]
% h2      : height of conductor [m]
% x       : separation distance [m]
% e_g     : permittivity of earth [F/m] (originaxly not used)
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Zg_self : Self earth impedance [Ohm/m]

m0=4*pi*1e-7;
Zg_self=zeros(con,con);

for x=1:1:con
   for y=1:1:con
       if x~=y
            D1=sqrt((h(1,x)-h(1,y))^2+d(x,y)^2);
            H=h(1,x)+h(1,y);
            p=1/sqrt(1i*omega*m0.*sigma_g);
            theta=rad2deg(atan(d(x,y)/H));
            if theta <= 50.45
                A=0.07360;
                a=0.1500;
            else
                A=0.00247*theta-0.05127;
                a=0.004726*theta-0.08852;
            end    

            % Self Impedance
            beta=(1-A*a)/(1-A);
            term1=A*log(sqrt((H+2*a*p).^2+d(x,y)^2)/D1);
            term2=(1-A)*log(sqrt((H+2*beta*p).^2+d(x,y)^2)/D1);

            Zg_self(x,y)=1i*omega.*m0/(2*pi).*(term1+term2);
       end
   end
end