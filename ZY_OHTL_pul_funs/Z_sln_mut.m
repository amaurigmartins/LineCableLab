% Z_snd_mut.m
function [Zg_mutual]=Z_sln_mut(h,d,e_g,m_g,sigma_g,omega,con)
% Function for the Mutual Earth Impedance by Semlyen's Approximate
% Formula

% Inputs
% h1      : height of conductor [m]
% h2      : height of conductor [m]
% d       : separation distance [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Zg_mutual: Mutual earth impedance [Ohm/m]

% Mutual Impedance

Zg_mutual=zeros(con,con);

for x=1:con
    for y=x+1:con
        if x~=y
            ter1ar_mut=(h(1,x)+h(1,y));
            ter1par_mut=pi*(d(x,y)^2+(h(1,x)+h(1,y))^2);
            ter1_mut=ter1ar_mut./ter1par_mut;

            ter2ar_mut=1i*omega.*m_g;
            ter2par_mut=sigma_g+1i*omega.*e_g;
            ter2_mut=sqrt(ter2ar_mut./ter2par_mut);

            Zg_mutual(x,y)=ter1_mut.*ter2_mut;
            Zg_mutual(y,x)=Zg_mutual(x,y);
        end
    end
end
