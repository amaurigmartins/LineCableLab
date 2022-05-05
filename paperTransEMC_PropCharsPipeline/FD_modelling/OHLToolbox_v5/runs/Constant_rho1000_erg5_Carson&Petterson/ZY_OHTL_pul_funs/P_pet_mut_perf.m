% Z_pet_slf.m
function [Ypg_mut]=P_pet_mut_perf(h,d,con)
% This function calculates the mutual potential coefficient (perfect ground)

% inputs:
% h1: height of conductor
% h2: height of conductor
% d12: distance of conductors
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Ppg_mutual

% Mutual potential coefficient
m0=4*pi*1e-7;
e0=8.854187817e-12;

Ypg_mut=zeros(con,con);


for x=1:1:con
    for y=1:1:con
        if x~=y
            % Mutual Impedance
            d2=sqrt(d(x,y).^2+(h(1,x)+h(1,y)).^2);
            d1=sqrt(d(x,y).^2+(h(1,x)-h(1,y)).^2);

            Ypg_mut(x,y)=log(d2./d1);
        end
    end
end
