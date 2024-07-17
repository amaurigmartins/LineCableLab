% Z_snd_mut.m
function [Zg_mutual]=Z_snd_mut(h,d,e_g,m_g,sigma_g,omega,con)
% Function for the Mutual Earth Impedance by Sunde's Approximate
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

m0=4*pi*1e-7;
gamma_g=sqrt(1i*omega.*m_g.*(sigma_g+1i*omega.*e_g));
Zg_mutual=zeros(con,con);
% Mutual Impedance
for x=1:con
    for y=x+1:con
      if x~=y
        ar_mut=(1+gamma_g.*((h(1,x)+h(1,y))/2)).^2+(gamma_g.*(d(x,y)/2)).^2;
        par_mut=(gamma_g.*((h(1,x)+h(1,y))/2)).^2+(gamma_g.*(d(x,y)/2)).^2;

        Zg_mutual(x,y)=1i*omega.*m0/(4*pi).*log(ar_mut./par_mut);
        Zg_mutual(y,x)=Zg_mutual(x,y);

      end
    end
end
