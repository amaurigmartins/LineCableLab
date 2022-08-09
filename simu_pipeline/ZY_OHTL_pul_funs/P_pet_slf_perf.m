% Z_pet_slf.m
function [Ypg_self]=P_pet_slf_perf(h,r,con)
% This function Calculates the self potential coefficient (perfect ground)

% inputs:
% h1: height of conductor
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Ppg_self

Ypg_self=zeros(con,con);

for k=1:1:con
    % Self potential coefficient
    Ypg_self(k,k)=log(2.*h(1,k)/r(k));
end