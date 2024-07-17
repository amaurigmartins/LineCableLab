% Z_pet_mut.m
function [Zg_mutual]=Z_pet_mut(h,d,e_g,m_g,sigma_g,omega,con)
% Function for the Mutual Earth Impedance by Pettersson's Approximate
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
b=sqrt(1i*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));

Zg_mutual=zeros(con,con);

% Mutual Impedance
for x=1:con
    for y=x+1:con
        if x~=y
            dp=sqrt((h(1,x)+h(1,y)+2./b).^2+d(x,y).^2);
            d2=sqrt(d(x,y).^2+(h(1,x)+h(1,y)).^2);
            term=log(dp./d2);
            Zg_mutual(x,y)=1i*omega.*m0./(2*pi).*term;
            Zg_mutual(y,x)=Zg_mutual(x,y);
        end
    end
end

