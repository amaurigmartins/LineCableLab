% Z_snd_mut.m
function [Zg_mutual]=Z_der_mut(h,d,m_g,sigma_g,omega,con)
% Function for the Self Earth Impedance by Alvarado - Betancourt Method
% Approximate Formula

% Inputs
% h1      : height of conductor [m]
% h2      : height of conductor [m]
% d       : separation distance [m]
% e_g     : permittivity of earth (originally not used) [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]
% Output
% Z_mutual: Mutual earth impedance [Ohm/m]

m0=4*pi*1e-7;
p=1./sqrt(1i*omega.*m_g*sigma_g);
%p=1./sqrt(j*omega.*m_g.*(sigma_g+j*omega.*e_g));
Zg_mutual=zeros(con,con);

if con == 1; return; end

% Mutual Impedance

for x=1:con
    for y=x+1:con
        if x==y
            Zg_mutual(x,y)=0;
        else

            terar_mut=sqrt((h(1,x)+h(1,y)+2*p).^2+d(x,y)^2);
            terpar_mut=sqrt((h(1,x)+h(1,y))^2+d(x,y)^2);
            Zg_mutual(x,y)=1i*omega.*m0./(2*pi).*log(terar_mut./terpar_mut);
            Zg_mutual(y,x)=Zg_mutual(x,y);
        end
    end
end
