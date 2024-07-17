% Z_snd_mut.m
function [Zg_mutual]=Z_alvadetan_mut(height,d,e_g,m_g,sigma_g,omega,con)
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
Zg_mutual=zeros(con,con);

for x=1:con
    for y=x+1:con
         if x~=y
             
            beta=d(x,y)./(height(1,x)+height(1,y));
            h=(height(1,x)+height(1,y))/2;

            ar11=(1+p./h).^2+beta.^2;
            par11=1+beta.^2;
            par12=(h./p.*(1+1i.*beta)+1).^3;
            par13=(h./p.*(1-1i.*beta)+1).^3;

            term1=0.25.*log(ar11./par11);
            term2=-1/48*(1./par12+1./par13);
            Jm=term1+term2;
            Zg_mutual(x,y)=1i*omega.*m0./(pi).*Jm;
            Zg_mutual(y,x)=Zg_mutual(x,y);

         end
    end
end
