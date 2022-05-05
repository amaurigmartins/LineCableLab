% Energization
function [vo1,Vo1]=src_enrg(amp,t_cl,data_t_sim,e)

f=50; % Frequency
o=amp*(cos(2*pi*f*data_t_sim+pi)); % Το συνημίτονο
samp=max(size(o)); % Size of samples

vo1=[zeros(1,int32(t_cl*e)),o(int32(t_cl*e+1):samp)]; % Το "αποκομμένο" συνημίτονο - int32: Signed 32-bit integer from -2,147,483,648 to 2,147,483,647

Vo1=fft(vo1); % FFT "αποκομμένου" συνημιτόνου (με το είδωλο)
Vo1=Vo1(1,1:(samp/2+1)); % FFT "αποκομμένου" συνημιτόνου (χωρίς το είδωλο)