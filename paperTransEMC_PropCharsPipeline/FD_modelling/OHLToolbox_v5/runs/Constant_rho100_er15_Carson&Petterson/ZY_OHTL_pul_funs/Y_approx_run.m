% Y_approx.m
% Simple Approximate Expressions
format long
close all
clear;
clc;


%% Input Data
% 1) Frequency
f=load('frequencies.txt'); % frequencies file
siz=length(f);

% 2) Overhead Transimission Line Configuration
h1=8.00;    % height of 1st conductor
h2=8.00;    % height of 2nd conductor
h3=8.00;    % height of 3rd conductor
d12=1.25;   % distance between 1st and 2nd conductor
d13=2.5;    % distance between 1st and 3rd conductor
d23=d12;    % distance between 2nd and 3rd conductor
r=0.015/2;  % external radius of conductor
rad_ex=r;
rad_in=0.0056/2; % internal radius of conductor

% 3) Electric Parameters
erg=10;         %relative permittivity of earth
mrg=1;          %relative permeability of earth
sigma_g=1/1000; %conductivity of earth
sigma_w=1/3.86063e-008; %conductivity of conductor

% Variables
e0=8.854187817e-12;
m0=4*pi*1e-7;
omega_total=2*pi*f;
e_g=e0*erg;
m_g=m0*mrg;

for k=1:siz
    omega=2*pi*f(k);
    n1_tetragwno=1;

%%%% Pettersson's Method %%%%

  %% Potential coefficients (self terms)

    %Potential Coeficients of Self Admittance - Perfect Ground
    Ps11_pet_perf=P_pet_slf_perf(h1,r);
    Ps22_pet_perf=P_pet_slf_perf(h2,r);
    Ps33_pet_perf=P_pet_slf_perf(h3,r);

    %Potential Coeficients of Self Admittance - Imperfect Ground
    Ps11_pet_imperf=P_pet_slf_imperf(h1,e_g,m_g,sigma_g,omega);
    Ps22_pet_imperf=P_pet_slf_imperf(h2,e_g,m_g,sigma_g,omega);
    Ps33_pet_imperf=P_pet_slf_imperf(h3,e_g,m_g,sigma_g,omega);
    
    % Total self
    Ps11_pet=Ps11_pet_perf+Ps11_pet_imperf;
    Ps22_pet=Ps22_pet_perf+Ps22_pet_imperf;
    Ps33_pet=Ps33_pet_perf+Ps33_pet_imperf;

    MN_self(k)=Ps11_pet_imperf;
    
  %% Potential coefficients (mutual terms)
    %Potential Coeficients of Mutal Admittance - Perfect Ground
    Pm12_pet_perf=P_pet_mut_perf(h1,h2,d12);
    Pm13_pet_perf=Y_pet_mut_perf(h1,h3,d13);
    Pm23_pet_perf=P_pet_mut_perf(h2,h3,d23);
    
    %Potential Coeficients of Mutal Admittance - Imperfect Ground
    [Pm12_pet_imperf,dq12]=P_pet_mut_imperf(h1,h2,d12,e_g,m_g,sigma_g,omega);
    [Pm13_pet_imperf,dq13]=Y_pet_mut_imperf(h1,h2,d13,e_g,m_g,sigma_g,omega);
    [Pm23_pet_imperf,dq23]=P_pet_mut_imperf(h2,h3,d23,e_g,m_g,sigma_g,omega);
  
    % Total Mutual
    Pm12_pet=Pm12_pet_perf+Pm12_pet_imperf;
    Pm13_pet=Pm13_pet_perf+Pm13_pet_imperf;
    Pm23_pet=Pm23_pet_perf+Pm23_pet_imperf;

    % Total Pottential Coefficient Matrix

    L_Q_mat=[Ps11_pet Pm12_pet Pm13_pet;Pm12_pet Ps22_pet Pm23_pet; Pm13_pet Pm23_pet Ps33_pet];

    % Total Admittance Matrix
    Ytot_Pet(:,:,k)=1i.*omega.*e0.*2.*pi.*n1_tetragwno.*inv(L_Q_mat);
    
end