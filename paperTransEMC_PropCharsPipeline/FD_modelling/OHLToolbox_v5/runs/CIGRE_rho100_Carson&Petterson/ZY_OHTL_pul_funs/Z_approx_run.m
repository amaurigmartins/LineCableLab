% Z_approx.m
% Simple Approximate Expressions
% format long
% close all
% clear;
% clc;

%% Input Data
% 1) Frequency
f=load('frequencies.txt'); % frequencies file
siz=length(f);

% 2) Overhead Transimission Line Configuration
h1=8.00;    % height of 1st conductor
h2=8.00;    % height of 2nd conductor
h3=8.00;    % height of 3rd conductor
d12=1.5;   % distance between 1st and 2nd conductor
d13=2.3;    % distance between 1st and 3rd conductor
d23=0.8;    % distance between 2nd and 3rd conductor
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

% Dynamic allocation 
Ztot_Carson=zeros(3,3,siz);
Ztot_Noda=zeros(3,3,siz);
Ztot_Deri=zeros(3,3,siz);
Ztot_Sunde=zeros(3,3,siz);
Ztot_Pettersson=zeros(3,3,siz);
Ztot_Semlyen=zeros(3,3,siz);

%% Series impedances
for k=1:siz
    omega=omega_total(k);
    freq=f(k);
%%%% Carson's Method %%%%
    % Formulas includes the effect of the perfect and the imperfect earth
    % Self Impedances
    Zs11_carson=Z_carson_slf(h1,rad_ex,freq,sigma_g,0);
    Zs22_carson=Z_carson_slf(h2,rad_ex,freq,sigma_g,0);
    Zs33_carson=Z_carson_slf(h3,rad_ex,freq,sigma_g,0);

    % Mutual Impedances
    Zm12_carson=Z_carson_mut(h1,h2,d12,freq,sigma_g,0);
    Zm13_carson=Z_carson_mut(h1,h3,d13,freq,sigma_g,0);
    Zm23_carson=Z_carson_mut(h2,h3,d23,freq,sigma_g,0);

%%%% Noda's Method %%%%
    % Formulas includes the effect of the perfect and the imperfect earth
    % Self Impedances
    Zs11_noda=Z_noda_slf(h1,0,sigma_g,omega,rad_ex);
    Zs22_noda=Z_noda_slf(h1,0,sigma_g,omega,rad_ex);
    Zs33_noda=Z_noda_slf(h1,0,sigma_g,omega,rad_ex);

    % Mutual Impedances
    Zm12_noda=Z_noda_mut(h1,h2,e_g,sigma_g,omega,d12);
    Zm13_noda=Z_noda_mut(h1,h3,e_g,sigma_g,omega,d13);
    Zm23_noda=Z_noda_mut(h2,h2,e_g,sigma_g,omega,d23);
    
%%%% Deri's Method %%%%
    
    % Self Impedances
    Zs11_deri=Z_der_slf(h1,e_g,m_g,sigma_g,omega);
    Zs22_deri=Z_der_slf(h2,e_g,m_g,sigma_g,omega);
    Zs33_deri=Z_der_slf(h3,e_g,m_g,sigma_g,omega);

    % Mutual Impedances
    Zm12_deri=Z_der_mut(h1,h2,d12,e_g,m_g,sigma_g,omega);
    Zm13_deri=Z_der_mut(h1,h3,d13,e_g,m_g,sigma_g,omega);
    Zm23_deri=Z_der_mut(h2,h3,d23,e_g,m_g,sigma_g,omega);
    
%%%% Sunde's Method %%%%
    
    % Self Impedances
    Zs11_sund=Z_snd_slf(h1,e_g,m_g,sigma_g,omega);
    Zs22_sund=Z_snd_slf(h2,e_g,m_g,sigma_g,omega);
    Zs33_sund=Z_snd_slf(h3,e_g,m_g,sigma_g,omega);

    % Mutual Impedances
    Zm12_sund=Z_snd_mut(h1,h2,d12,e_g,m_g,sigma_g,omega);
    Zm13_sund=Z_snd_mut(h1,h3,d13,e_g,m_g,sigma_g,omega);
    Zm23_sund=Z_snd_mut(h2,h3,d23,e_g,m_g,sigma_g,omega);


%%%% Semlyen's Method %%%%
    
    % Self Impedances
    Zs11_seml=Z_sln_slf(h1,e_g,m_g,sigma_g,omega);
    Zs22_seml=Z_sln_slf(h2,e_g,m_g,sigma_g,omega);
    Zs33_seml=Z_sln_slf(h3,e_g,m_g,sigma_g,omega);

    % Mutual Impedances
    Zm12_seml=Z_sln_mut(h1,h2,d12,e_g,m_g,sigma_g,omega);
    Zm13_seml=Z_sln_mut(h1,h3,d13,e_g,m_g,sigma_g,omega);
    Zm23_seml=Z_sln_mut(h2,h3,d23,e_g,m_g,sigma_g,omega);


%%%% Pettersson's Method %%%%

    % Self Impedances
    Zs11_pet=Z_pet_slf(h1,e_g,m_g,sigma_g,omega);
    Zs22_pet=Z_pet_slf(h2,e_g,m_g,sigma_g,omega);
    Zs33_pet=Z_pet_slf(h3,e_g,m_g,sigma_g,omega);

    % Mutual Impedances
    Zm12_pet=Z_pet_mut(h1,h2,d12,e_g,m_g,sigma_g,omega);
    Zm13_pet=Z_pet_mut(h1,h3,d13,e_g,m_g,sigma_g,omega);
    Zm23_pet=Z_pet_mut(h2,h3,d23,e_g,m_g,sigma_g,omega);

%%%% Total Matrices %%%%
    Zg_Carson=[Zs11_carson Zm12_carson Zm13_carson; Zm12_carson Zs22_carson Zm23_carson; Zm13_carson Zm23_carson Zs33_carson];
    Zg_Noda=[Zs11_noda Zm12_noda Zm13_noda; Zm12_noda Zs22_noda Zm23_noda; Zm13_noda Zm23_noda Zs33_noda];
    Zg_Deri=[Zs11_deri Zm12_deri Zm13_deri; Zm12_deri Zs22_deri Zm23_deri; Zm13_deri Zm23_deri Zs33_deri];
    Zg_Sund=[Zs11_sund Zm12_sund Zm13_sund; Zm12_sund Zs22_sund Zm23_sund; Zm13_sund Zm23_sund Zs33_sund];
    Zg_Pet=[Zs11_pet Zm12_pet Zm13_pet; Zm12_pet Zs22_pet Zm23_pet; Zm13_pet Zm23_pet Zs33_pet];
    Zg_Sem=[Zs11_seml Zm12_seml Zm13_seml; Zm12_seml Zs22_seml Zm23_seml; Zm13_seml Zm23_seml Zs33_seml];
       
    Zpg=Z_self_mut_pg(h1,h2,h3,d12,d23,d13,r,omega);        % Influence of perfect earth
    Zskin_self=skeffct_tb_fun(rad_ex,rad_in,sigma_w,omega); % Skin effect
    
    Zskin=[Zskin_self 0 0; 0 Zskin_self 0; 0 0 Zskin_self];
    
    Ztot_Carson(:,:,k)=Zskin+Zg_Carson;
    Ztot_Noda(:,:,k)=Zskin+Zg_Noda;
    Ztot_Deri(:,:,k)=Zskin+Zpg+Zg_Deri;
    Ztot_Sunde(:,:,k)=Zskin+Zpg+Zg_Sund;
    Ztot_Pettersson(:,:,k)=Zskin+Zpg+Zg_Pet;
    Ztot_Semlyen(:,:,k)=Zskin+Zpg+Zg_Sem;    
end

%% Plot results
% Self impedance (Z11)
figure(1)
subplot(2,1,1)
loglog(f,squeeze(abs(Ztot_Carson(1,1,:))),f,squeeze(abs(Ztot_Noda(1,1,:))),f,squeeze(abs(Ztot_Deri(1,1,:))),f,squeeze(abs(Ztot_Sunde(1,1,:))),f,squeeze(abs(Ztot_Pettersson(1,1,:))),f,squeeze(abs(Ztot_Semlyen(1,1,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Carson','Noda','Deri','Sunde','Pettersson','Semlyen')
grid
title('Self impedance - Z11')

subplot(2,1,2)
loglog(f,squeeze(angle(Ztot_Carson(1,1,:))),f,squeeze(angle(Ztot_Noda(1,1,:))),f,squeeze(angle(Ztot_Deri(1,1,:))),f,squeeze(angle(Ztot_Sunde(1,1,:))),f,squeeze(angle(Ztot_Pettersson(1,1,:))),f,squeeze(angle(Ztot_Semlyen(1,1,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Carson','Noda','Deri','Sunde','Pettersson','Semlyen')
grid

% Mutual impedance (Z12)
figure(2)
subplot(2,1,1)
loglog(f,squeeze(abs(Ztot_Carson(1,2,:))),f,squeeze(abs(Ztot_Noda(1,2,:))),f,squeeze(abs(Ztot_Deri(1,2,:))),f,squeeze(abs(Ztot_Sunde(1,2,:))),f,squeeze(abs(Ztot_Pettersson(1,2,:))),f,squeeze(abs(Ztot_Semlyen(1,2,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Carson','Noda','Deri','Sunde','Pettersson','Semlyen')
grid
title('Mutual impedance - Z12')

subplot(2,1,2)
loglog(f,squeeze(angle(Ztot_Carson(1,2,:))),f,squeeze(angle(Ztot_Noda(1,2,:))),f,squeeze(angle(Ztot_Deri(1,2,:))),f,squeeze(angle(Ztot_Sunde(1,2,:))),f,squeeze(angle(Ztot_Pettersson(1,2,:))),f,squeeze(angle(Ztot_Semlyen(1,2,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Carson','Noda','Deri','Sunde','Pettersson','Semlyen')
grid

% Mutual impedance (Z13)
figure(3)
subplot(2,1,1)
loglog(f,squeeze(abs(Ztot_Carson(1,3,:))),f,squeeze(abs(Ztot_Noda(1,3,:))),f,squeeze(abs(Ztot_Deri(1,3,:))),f,squeeze(abs(Ztot_Sunde(1,3,:))),f,squeeze(abs(Ztot_Pettersson(1,3,:))),f,squeeze(abs(Ztot_Semlyen(1,3,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Carson','Noda','Deri','Sunde','Pettersson','Semlyen')
grid
title('Mutual impedance - Z13')

subplot(2,1,2)
loglog(f,squeeze(angle(Ztot_Carson(1,3,:))),f,squeeze(angle(Ztot_Noda(1,3,:))),f,squeeze(angle(Ztot_Deri(1,3,:))),f,squeeze(angle(Ztot_Sunde(1,3,:))),f,squeeze(angle(Ztot_Pettersson(1,3,:))),f,squeeze(angle(Ztot_Semlyen(1,3,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Carson','Noda','Deri','Sunde','Pettersson','Semlyen')
grid

