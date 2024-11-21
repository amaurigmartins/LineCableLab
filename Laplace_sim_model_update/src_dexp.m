% Double exponential
function [vo1,Vo1,c]=src_dexp(amp,time_cl_brkr,time_sim,samples,e,f,data_t_sim)

t_surge=time_sim-time_cl_brkr; % Time after the switch closing
data_srg=0:1/e:t_surge; % Vector of time for the double exponential

%c=0; % FFT
%c=2*2*pi*(1/time_sim); % Wilcox
c=log(samples^2)/time_sim; % Wedepohl

%s=c+1i*2*pi*f;

%alpha=72999.79; % 0.1/10 ìs
%beta=68516774.89;

%alpha=36499.9; % 0.1/20 ìs
%beta=76515953.19;

%alpha=24333.26; % 0.1/30 ìs
%beta=81145936.42;

%alpha=18249.95; % 0.1/40 ìs
%beta=84411179.69;

%alpha=14600; % 0.1/50 ìs
%beta=86933087.6;

%alpha=14600; % 0.2/50 ìs
%beta=39534073.548;

%alpha=14600; % 0.4/50 ìs
%beta=17776060.531;

%alpha=14600; % 0.6/50 ìs
%beta=11065546.633;

%alpha=14600; % 0.8/50 ìs
%beta=7878127.262;

%alpha=14600; % 1/50 ìs
%beta=6039683.134;

%alpha=288000; % 1.2/5 ìs
%beta=1241666.667;

alpha=14600; % 1.2/50 ìs
beta=2466666.667;
corr_fact=2.5/2.41;
%alpha=3500; % 1.2/200 ìs
%beta=2625000;

%alpha=730; % 1.2/1000 ìs
%beta=7722951.77;

%alpha=292; % 1.2/2500 ìs
%beta=68573118.46;

%alpha=291.9991617; % 250/2500 ìs
%beta=16406.80018;

d_emp_fun=amp*corr_fact*((exp(-alpha*data_srg))-(exp(-beta*data_srg)));

vo1=[zeros(1,int32(time_cl_brkr*e)),d_emp_fun]; % Double exponential with switch - int32: Signed 32-bit integer from -2,147,483,648 to 2,147,483,647

Cn=exp(-c*data_t_sim)*(1/e);

Vo1=fft(vo1.*Cn); % FFT of double exponential (with image)
Vo1=Vo1(1,1:(round(samples/2)+1)); % FFT of double exponential (without image)
