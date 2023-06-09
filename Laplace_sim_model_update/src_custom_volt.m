function [vo1,Vo1,c]=src_custom_volt(t_sim,samples,e,data_t_sim)

%c=0; % FFT
%c=2*2*pi*(1/time_sim); % Wilcox
c=log(samples^2)/t_sim; % Wedepohl

load source_A50.mat

time_max_field=t(length(t));
rest_t=t_sim-time_max_field;
data_t_surge=0:1/e:time_max_field; % Vector of time

Smooth_Csaps=csaps(t,Isrc);
Source_Csaps=ppval(Smooth_Csaps,data_t_surge);

vo1=[Source_Csaps,zeros(1,int32(rest_t*e))];
% vo1=vo1/max(vo1);

Cn=exp(-c*data_t_sim)*(1/e);

Vo1=fft(vo1.*Cn); % FFT double exponential (with image)
Vo1=Vo1(1,1:(samples/2+1)); % FFT double exponential (without image)