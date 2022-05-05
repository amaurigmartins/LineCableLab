% P_wise_slf_imperf.m
function [Pg_self]=P_wise_slf_imperf(h,e_g,m_g,sigma_g,omega,con)
% This function Calculates the self potential coefficient (imperfect ground - Wise's Integral formula)

% inputs:
% h1: height of conductor
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Pg_self

m0=4*pi*1e-7;
e0=8.854187817e-12;
k0=omega*(sqrt(m0*e0));
gama0 = 1i*omega*sqrt(m0*e0);
gama_g = sqrt(1i*omega*m0*(sigma_g + 1i*omega*e_g));

Pg_self=zeros(con,con);

for k=1:1:con

    % Self potential coefficient
    H=2*h(1,k);

    y_lambda =@(lambda) sum([0 (((2.*m_g.*(gama0.^2).*(m0.*lambda + (sqrt(lambda.^2 + gama_g.^2 + k0.^2)).*m_g).*exp(-H.*lambda))/(((sqrt(lambda.^2 + gama_g.^2 + k0.^2)).*m0 + lambda.*m_g).*((sqrt(lambda.^2 + gama_g.^2 + k0.^2)).*(gama0.^2).*m_g + lambda.*(gama_g.^2).*m0))))], 'omitnan');
    %     y_lambda=@(lambda) sum([0 y_lambda(lambda)],'omitnan');
    Qm = integral(y_lambda,0,Inf,'ArrayValued',true);
%     if isinf(abs(Qm))
%         Qm = complex(1/eps,1/eps);
%     end
    Pg_self(k,k)=Qm;
end