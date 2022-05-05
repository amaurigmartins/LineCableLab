% Z_pet_slf.m
function [Ypg_mut]=Y_pet_mut_perf(h1,h2,d12)
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

% Mutual Impedance
d2=sqrt(d12.^2+(h1+h2).^2);
d1=sqrt(d12.^2+(h1-h2).^2);

Ypg_mut=log(d2./d1);
