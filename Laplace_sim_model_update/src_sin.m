% Sinusoidal Voltage
function [vo1,Vo1]=src_sin(amp,data_t_sim)

f=50; % Frequency
vo1=amp*(cos(2*pi*f*data_t_sim)); % cosine
samp=max(size(data_t_sim)); % Size of samples

Vo1=fft(vo1); % FFT cosine (with image)
Vo1=Vo1(1,1:(samp/2+1)); % FFT cosine (without image)