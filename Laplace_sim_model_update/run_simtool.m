%% Simulation Tool (SimTool)
close all; % Closes all figures whose handles are not hidden
clear; % Remove items from workspace, freeing up system memory
clc; % Clears all input and output from the Command Window display, giving you a "clean screen"
format longEng; % Engineering format that has exactly 16 significant digits and a power that is a multiple of three.

%% Input Data (1)
ord=6; % Number of conductors
file='ZY_Carson_CP100.mat'; % pul parameter input file
time_sim=0.1; % Simulation time (sec)
time_cl_brkr=0; % Time close of breaker (sec)
length_line=1000; % Line length (m)

%% Simulation Options (2)
sim_flag=4; % Frequency scan = 0, TD simulation >=1
% a) Sinusoidal Voltage, flag=1
% b) Energization, flag=2
% c) Step Response, flag=3
% d) Double Exponential Voltage Source, flag=4
% e) Custom Source Voltage, flag=5
% f) 3ph Source, flag=6
% g1) LI Current Source, flag=71
% g2) A50 Current Source, flag=72
% g3) F5 Current Source, flag=73
% g4) F95 Current Source, flag=74
% g5) S5 Current Source, flag=75
% g6) S95 Current Source, flag=76

amp=1; % Source amplitude in pu!!

e=1.5*5000000; % Sampling rate. Default: e=10000000;
data_t_sim=0:1/e:time_sim; % Vector of time
samples=max(size(data_t_sim));  % Size of samples in time domain
f=0:e/samples:(e/2); % Vector of frequency

%%  ZY parameters files (3)
%[gamma_dis,Z_dis,Y_dis,Ti_dis,freq]=pul_clc(ord,num_files); % Function pul_clc - Calculates matrix Y(num_files x ord), matrices Z' and Y' (num_files x ord^2) and matrix Ti (num_files x ord^2) at each freq of the total num_files
[gamma_dis,Z_dis,Y_dis,Ti_dis,freq,num]=pul_input_fun(file);
num_files=num; % Number of input pul data
%[Zmodal_ch,gamma,Z,Y,Ti,freq]=pul_clc_finite_length(ord,num_files,length);

%% Build Circuit Struct - Termination Conditions (4A)
Ys_dis=zeros(num_files,ord);
Yr_dis=zeros(num_files,ord);

for o=1:1:num_files
% Set up simulations - termination
   Ys_dis(o,:)=[0+1j*2*pi*freq(o)*0 (1/574)+1j*2*pi*freq(o)*0 (1/574)+1j*2*pi*freq(o)*0 (1/306)+1j*2*pi*freq(o)*0 (1/1354)+1j*2*pi*freq(o)*0 (1/1354)+1j*2*pi*freq(o)*0];
   Yr_dis(o,:)=[(1/574)+1j*2*pi*freq(o)*0 (1/574)+1j*2*pi*freq(o)*0 (1/574)+1j*2*pi*freq(o)*0 (1/306)+1j*2*pi*freq(o)*0 (1/1354)+1j*2*pi*freq(o)*0 (1/1354)+1j*2*pi*freq(o)*0];
end
%% Simulation (5)
if sim_flag >= 1 % A)Time-Domain Simulation
    [v]=td_sim(ord,f,freq,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis,length_line,time_cl_brkr,samples,e,data_t_sim,time_sim,amp,sim_flag);
else % B)Frequency scan
    Vo1=freq_scan_sim(ord,freq,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis,length_line);
end

%% Save Results (6)
%save fs_WiseCG_len100.mat Vo1 freq
save LI_Carson_CP100_1000.mat v data_t_sim