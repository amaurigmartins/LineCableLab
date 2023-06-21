% Energization
function [vo1,Vo1,c]=src_enrg(amp,t_cl,data_t_sim,e)

c=0; % for FFT
f=50; % Frequency
o=amp*(cos(2*pi*f*data_t_sim+pi)); % cosine
samp=max(size(o)); % Size of samples

vo1=[zeros(1,int32(t_cl*e)),o(int32(t_cl*e+1):samp)]; % "Clipped" cosine - int32: Signed 32-bit integer from -2,147,483,648 to 2,147,483,647

Vo1=fft(vo1); % FFT of the "clipped" cosine (with image)
Vo1=Vo1(1,1:(samp/2+1)); % FFT of the "clipped" cosine (without image)