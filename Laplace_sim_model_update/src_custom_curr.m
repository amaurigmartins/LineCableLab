function [io1,Io1,c]=src_custom_curr(t_sim,samples,e,data_t_sim,sim_flag)

%% Options
%c=0; % FFT
%c=2*2*pi*(1/time_sim); % Wilcox
c=log(samples^2)/t_sim; % Wedepohl


%% Load sources
if sim_flag == 71 % LI current source, flag=71
    load('source_LI.mat');
elseif sim_flag == 72 % A50 current source, flag=72
    load('source_A50.mat');
elseif sim_flag == 73 % F5 current source, flag=73
    load('source_F5.mat');    
elseif sim_flag == 74 % F95 current source, flag=74
    load('source_F95.mat');    
elseif sim_flag == 75 % S5 current source, flag=75
    load('source_S5.mat');    
elseif sim_flag == 76 % S95 current source, flag=76
    load('source_S95.mat');
end

%% Generate FD Current Source
time_max_field=t(length(t));
rest_t=t_sim-time_max_field;
data_t_surge=0:1/e:time_max_field; % Vector of time

Smooth_Csaps=csaps(t,Isrc);
Source_Csaps=ppval(Smooth_Csaps,data_t_surge);

io1=[Source_Csaps,zeros(1,int32(rest_t*e))];
% vo1=vo1/max(vo1);

Cn=exp(-c*data_t_sim)*(1/e);

Io1=fft(io1.*Cn); % FFT with image
Io1=Io1(1,1:(samples/2+1)); % FFT without image