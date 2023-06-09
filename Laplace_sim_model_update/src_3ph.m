% 3ph
function [vo1,vo2,vo3,Vo1,Vo2,Vo3,c]=src_3ph(amp,t_cl,data_t_sim,e,samples,f,time_sim)

c=0; % FFT
%c=2*2*pi*(1/time_sim); % Wilcox
%c=log(samples^2)/time_sim; % Wedepohl

%s=c+1i*2*pi*f;

f=50; % Frequency
o1=amp*(cos(2*pi*f*data_t_sim)); % Cosine of phase a
o2=amp*(cos(2*pi*f*data_t_sim-(2*pi/3))); % Cosine of phase b
o3=amp*(cos(2*pi*f*data_t_sim+(2*pi/3))); % Cosine of phase c

vo1=[zeros(1,int32(t_cl*e)),o1(int32(t_cl*e+1):samples)]; % "Choked" cosine of phase a - int32: Signed 32-bit integer from -2,147,483,648 to 2,147,483,647
vo2=o2;
vo3=o3;

Cn=exp(-c*data_t_sim)*(1/e);

Vo1=fft(vo1.*Cn); % FFT of the "choked" cosine (with image)
Vo2=fft(vo2.*Cn);
Vo3=fft(vo3.*Cn);

Vo1=Vo1(1,1:(samples/2+1)); % FFT of the "choked" cosine (without image)
Vo2=Vo2(1,1:(samples/2+1));
Vo3=Vo3(1,1:(samples/2+1));