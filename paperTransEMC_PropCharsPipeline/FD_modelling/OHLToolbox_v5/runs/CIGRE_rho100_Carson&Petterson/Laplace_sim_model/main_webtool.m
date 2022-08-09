%% Main Function
% Figs 1: vo1(t)
% Fig  2: vo2(t)
% Fig  3: v1(t)
% Fig  4: v2(t)
% ??? ????? ???? ?? ???????? ??????? ?? ord ??? ??????????...
close all; % Deletes all figures whose handles are not hidden
clear; % Remove items from workspace, freeing up system memory
clc; % Clears all input and output from the Command Window display, giving you a "clean screen"
format longEng; % Engineering format that has exactly 16 significant digits and a power that is a multiple of three.

%%
% 0) Data
ord=3; % Number of conductors
file='ZY_sim_ZY_Carson_CP1000.mat'; % pul parameter input file
time_sim=0.1; % Simulation time (sec)
time_cl_brkr=0; % Time close of breaker (sec)
length=100; % Line length (m)
sim_flag=0; % TD simulation = 0, Frequency scan = 1

%%
% 1) Sampling Rate
e=100000000; % Sampling rate e=10000000;
data_t_sim=0:1/e:time_sim; % Vector of time
samples=max(size(data_t_sim));  % Size of samples in time domain
f=0:e/samples:(e/2); % Vector of frequency

%%
% 2) pul Calculation
%[gamma_dis,Z_dis,Y_dis,Ti_dis,freq]=pul_clc(ord,num_files); % Function pul_clc - ?????????? ?? ?????? ? (num_files x ord), ???? ??????? ?' ??? Y' (num_files x ord^2) ??? ?? ?????? Ti (num_files x ord^2) ?? ???? ????????? freq ??? ??? num_files ??? ???????
[gamma_dis,Z_dis,Y_dis,Ti_dis,freq,num]=pul_input_fun(file);
num_files=num; % Number of input pul data
%[Zmodal_ch,gamma,Z,Y,Ti,freq]=pul_clc_finite_length(ord,num_files,length);

%% 3A) Circuit Struct - Termination Conditions
Ys_dis=zeros(num_files,ord);
Yr_dis=zeros(num_files,ord);

for o=1:1:num_files
%     Ys_dis(o,:)=[0+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0]; % O ??????? ???????????? ?? Siemens ?? ???? ???????? ????????? ??? ??? num_files ??? ?? ?????? ??? ???? S - (1 x ord)
%     Yr_dis(o,:)=[(1/500)+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0]; % O ??????? ???????????? ?? Siemens ?? ???? ???????? ????????? ??? ??? num_files ??? ?? ?????? ??? ???? R - (1 x ord)
    Ys_dis(o,:)=[0+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0];
    Yr_dis(o,:)=[(1/500)+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0 (1/500)+1j*2*pi*freq(o)*0];
end


%% 4) Simulation
if sim_flag == 0
    % A)Time-Domain Simulation
    v=td_sim(ord,f,freq,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis,length,time_cl_brkr,samples,e,data_t_sim,time_sim);
else
    % B)Frequency scan
    %[Vo1]=freq_scan(f); % Function src_sin
    Vo1=freq_scan_sim(ord,gamma_dis,Z_dis,Y_dis,Ti_dis,freq,Ys_dis,Yr_dis,length);
end
save FS_Noda_CP1000.mat Vo1 freq
%save Step_Carson_CP1000_100.mat v data_t_sim