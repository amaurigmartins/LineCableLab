%% Init
clear;
clc;
close all;
addpath(genpath('ohltoolbox'));

%% Frequency range  
points_per_dec=10;
f_dec1=1:10/points_per_dec:10;
f_dec2=10:100/points_per_dec:100;
% f_dec3=100:1000/points_per_dec:1000;
% f_dec4=1000:10000/points_per_dec:10000;
% f_dec5=10000:100000/points_per_dec:100000;
% f_dec6=100000:1000000/points_per_dec:1000000;
% f=transpose([f_dec1(1:length(f_dec1)-1) f_dec2(1:length(f_dec2)-1) f_dec3(1:length(f_dec3)-1) f_dec4(1:length(f_dec4)-1) f_dec5(1:length(f_dec5)-1) f_dec6]);

f=transpose([f_dec1(1:length(f_dec1)-1) f_dec2]);
freq_siz=length(f);

%% Line parameters
[line_length,ord,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d]=set_line_data();

%% Flags
ZYprnt=0; % Flag to print parameters
FD_flag=0; % Flag to calculate parameters using the Longmire-Smith FD approach
decmp_flag = 2; % Modal decomposition flag. (1 QR,(2)QR ATP,(3)simple_QR,(4)simple_QR,(5)NR,(6)NR_back,(7)SQP,(8)SQP_back,(9)LM,(10)LM_back,(11)LM_fast,(12)LM_alt

%% Calculations
tic
[Ztot_Carson,~,~,~,~,~,~] = Z_clc_fun(f,ord,ZYprnt,FD_flag,freq_siz,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d); % Calculate Z pul parameters by different earth approaches
Z=Ztot_Carson;
[Ytot_Imag,~] = Y_clc_fun(f,ord,ZYprnt,FD_flag,freq_siz,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d); % Calculate Y pul parameters by different earth approaches
Y=Ytot_Imag;
[Zch_mod,Ych_mod,Zch,Ych,g_dis,Ti_dis,Z_dis,Y_dis] = mode_decomp_fun(Z,Y,f,freq_siz,ord,decmp_flag); % Modal decomposition
[H_mod,F_mod,pol_co] = HF_VF_fun(Ti_dis,g_dis,line_length,f,ord); % Vector Fitting
toc

% syms s
% num=(pol_co(3,3).b);
% den=(pol_co(3,3).a);
% G = poly2sym(num,s)/poly2sym(den,s)
% G = matlabFunction(G)'
% 
% figure
% plot(f,G(f*2*pi))