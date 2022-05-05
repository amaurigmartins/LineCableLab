function [vo1,Vo1]=src_custom(t_sim,samp,e)

%load time_CS.mat
%load voltage_CS.mat

load time_SS.mat
load voltage_SS.mat

time_max_field=time_SS(length(time_SS));
rest_t=t_sim-time_max_field;
data_t_surge=0:1/e:time_max_field; % Vector of time

Smooth_Csaps=csaps(time_SS,voltage_SS);
Source_Csaps=ppval(Smooth_Csaps,data_t_surge);

vo1=[Source_Csaps,zeros(1,int32(rest_t*e)+1)];
vo1=vo1/max(vo1);

Vo1=fft(vo1); % FFT διπλεκθετικής (με το είδωλο)
Vo1=Vo1(1,1:(samp/2+1)); % FFT διπλεκθετικής (χωρίς το είδωλο)