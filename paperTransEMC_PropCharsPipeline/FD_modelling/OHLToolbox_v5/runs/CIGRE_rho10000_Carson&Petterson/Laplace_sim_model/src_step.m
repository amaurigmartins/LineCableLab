% Step response
function [vo1,Vo1,c]=src_step(amp,time_cl_brkr,samples,e,f,data_t_sim,time_sim)

%c=0; % FFT
c=2*pi*(1/time_sim); % Wilcox - ’λλαξε η τιμή για τον regular!!
%c=log((samples^2)/4)/time_sim; % Wedepohl - ’λλαξε η τιμή για τον regular!!

%s=c+1i*2*pi*f;

o=amp*ones(1,samples); % Η μοναδιαία "παύλα"
%vo1=[zeros(1,int32(time_cl_brkr*e)),o(int32(time_cl_brkr*e+1):samples)]; % Η step - int32: Signed 32-bit integer from -2,147,483,648 to 2,147,483,647
vo1=[zeros(1,int32(time_cl_brkr*e)),o(int32(time_cl_brkr*e+1):int32((time_sim-time_cl_brkr)*e)),zeros(1,int32(time_cl_brkr*e+1))]; % Η step - int32: Signed 32-bit integer from -2,147,483,648 to 2,147,483,647


Cn=exp(-c*data_t_sim)*(1/e);

Vo1=fft(vo1.*Cn); % FFT step (με το είδωλο)
Vo1=Vo1(1,1:(samples/2+1)); % FFT step (χωρίς το είδωλο)