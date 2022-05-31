%% Initialization
close all; % Deletes all figures whose handles are not hidden
clear; % Remove items from workspace, freeing up system memory
clc; % Clears all input and output from the Command Window display, giving you a "clean screen"
format longEng; % Engineering format that has exactly 16 significant digits and a power that is a multiple of three.

%% Library functions
currmfile = mfilename('fullpath');
currPath = currmfile(1:end-length(mfilename()));
addpath([currPath 'ZY_OHTL_pul_funs']);
addpath([currPath 'mode_decomp_funs']);
addpath([currPath 'LongmireSmith_funs']);

%% Frequency range 
points_per_dec=10;
f_dec1=1:10/points_per_dec:10;
f_dec2=10:100/points_per_dec:100;
f_dec3=100:1000/points_per_dec:1000;
f_dec4=1000:10000/points_per_dec:10000;
f_dec5=10000:100000/points_per_dec:100000;
f_dec6=100000:1000000/points_per_dec:1000000;
%f=transpose([f_dec1(1:length(f_dec1)-1) f_dec2(1:length(f_dec2)-1) f_dec3(1:length(f_dec3)-1) f_dec4(1:length(f_dec4)-1) f_dec5(1:length(f_dec5)-1) f_dec6]);
f=transpose([f_dec4(1:length(f_dec4)-1) f_dec5(1:length(f_dec5)-1) f_dec6]); % for NB-PLC
%f_dec7=1000000:10000:100000000;
%f=transpose(f_dec7); % for BB-PLC
%f=transpose([1E-6 f_dec1(1:length(f_dec1)-1) f_dec2(1:length(f_dec2)-1) f_dec3(1:length(f_dec3)-1) f_dec4(1:length(f_dec4)-1) f_dec5(1:length(f_dec5)-1) f_dec6(1:length(f_dec6)-1) f_dec7]);

freq_siz=length(f);

%% Line Parameters
[line_length,ord,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d]=LineData_fun();

%% Flags
ZYprnt=1; % Flag to print parameters
FD_flag=0; % Flag to calculate parameters using the Longmire-Smith FD approach
decmp_flag = 5; % Modal decomposition flag. (1 QR,(2)QR ATP,(3)simple_QR,(4)simple_QR,(5)NR,(6)NR_back,(7)SQP,(8)SQP_back,(9)LM,(10)LM_back,(11)LM_fast,(12)LM_alt

%% Calculations
tic
[Ztot_Carson,Ztot_Noda,Ztot_Deri,Ztot_AlDe,Ztot_Sunde,Ztot_Pettersson,Ztot_Semlyen] = Z_clc_fun(f,ord,ZYprnt,FD_flag,freq_siz,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d); % Calculate Z pul parameters by different earth approaches
[Ytot_Imag,Ytot_Pettersson] = Y_clc_fun(f,ord,ZYprnt,FD_flag,freq_siz,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d); % Calculate Y pul parameters by different earth approaches
[Zch_mod,Ych_mod,Zch,Ych,g_dis,Ti_dis,Z_dis,Y_dis] = mode_decomp_fun(Ztot_AlDe,Ytot_Imag,f,freq_siz,ord,decmp_flag); % Modal decomposition
[H_mod,F_mod,pol_co] = HF_VF_fun(Ti_dis,g_dis,line_length,f,ord); % Vector Fitting
toc

%% Save files
%save ZY_LM_CP_BB_1000.mat Ztot_Sunde Ytot_Imag g_dis Zch_mod Z_dis Y_dis Ti_dis f
%save ZY_sim_ZY_AlDe_CP1000 g_dis Z_dis Y_dis Ti_dis f
