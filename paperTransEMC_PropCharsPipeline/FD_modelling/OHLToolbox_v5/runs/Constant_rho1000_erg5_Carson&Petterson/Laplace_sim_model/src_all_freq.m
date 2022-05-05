% Sinusoidal Voltage
function [Vo1,c]=src_all_freq(data_t_sim,samples,time_sim,f,e)

%c=0; % FFT
%c=2*2*pi*(1/time_sim); % Wilcox
c=log(samples^2)/time_sim; % Wedepohl

%s=c+1i*2*pi*f;

Cn=exp(-c*data_t_sim)*(1/e);

Vo1=ones(1,samples/2+1);