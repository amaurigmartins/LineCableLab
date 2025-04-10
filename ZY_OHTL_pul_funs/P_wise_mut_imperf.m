% P_wise_mut_imperf.m
function [Pg_mutual]=P_wise_mut_imperf(h,d,e_g,m_g,sigma_g,omega,con)
% This function Calculates the self potential coefficient (imperfect ground -  Wise's Integral formula)

% inputs:
% h1: height of conductor
% h2: height of conductor
% d12: distance of conductors
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Pg_mutual


% Mutual Impedance
m0=4*pi*1e-7;
e0=8.854187817e-12;
k0=omega*(sqrt(m0*e0));
gama0 = 1i*omega*sqrt(m0*e0);
gama_g = sqrt(1i*omega*m0*(sigma_g + 1i*omega*e_g));

Pg_mutual=zeros(con,con);

if con == 1; return; end


for x=1:con
    for y=x+1:con
        if x~=y
            
            h1=h(1,x);
            h2=h(1,y);
            H=h1+h2;
            d12 = d(x,y);
            
            y_lambda =@(lambda) sum([0 (((2*m_g*(gama0^2)*(m0*lambda + (sqrt(lambda^2 + gama_g^2 + k0^2))*m_g)*exp(-H*lambda))/(((sqrt(lambda^2 + gama_g^2 + k0^2))*m0 + lambda*m_g)*((sqrt(lambda^2 + gama_g^2 + k0^2))*(gama0^2)*m_g + lambda*(gama_g^2)*m0)))*cos(d12*lambda))],'omitnan');
%             y_lambda=@(lambda) sum([0 y_lambda(lambda)],'omitnan');
            Qm = integral(y_lambda,0,Inf,'ArrayValued',true);
            
            Pg_mutual(x,y)=Qm;
            Pg_mutual(y,x)=Pg_mutual(x,y);
        end
    end
end
