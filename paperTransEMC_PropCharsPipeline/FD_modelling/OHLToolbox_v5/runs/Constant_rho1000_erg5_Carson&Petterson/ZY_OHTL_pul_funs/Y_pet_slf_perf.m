% Z_pet_slf.m
function [Ypg_self]=Y_pet_slf_perf(h1,r)
% This function Calculates the self potential coefficient (perfect ground)

% inputs:
% h1: height of conductor
% e_g: permittivity of earth
% m_g: permeability of earth
% sigma_g: conductivity of earth
% omega: angular frequency

% output: Ppg_self

% Self potential coefficient
Ypg_self=log(2.*h1/r);
