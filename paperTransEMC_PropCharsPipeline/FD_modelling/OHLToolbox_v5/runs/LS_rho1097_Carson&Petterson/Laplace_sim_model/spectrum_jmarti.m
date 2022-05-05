e=20000000;
samples=max(size(t));
f_atp=0:e/samples:(e/2);

VRa=fft(vRa);
VRa=VRa(1:(samples/2+1),1);
VRa=abs(VRa)/(samples/2);

VRc=fft(vRc);
VRc=VRc(1:(samples/2+1),1);
VRc=abs(VRc)/(samples/2);

VM1a=fft(vM1a);
VM1a=VM1a(1:(samples/2+1),1);
VM1a=abs(VM1a)/(samples/2);

VM2b=fft(vM2b);
VM2b=VM2b(1:(samples/2+1),1);
VM2b=abs(VM2b)/(samples/2);
% 
% Vs1c=fft(vVs1c);
% Vs1c=Vs1c(1:(samples/2+1),1);
% Vs1c=abs(Vs1c)/(samples/2);
% 
% Vs2a=fft(vVs2a);
% Vs2a=Vs2a(1:(samples/2+1),1);
% Vs2a=abs(Vs2a)/(samples/2);
% 
% Vs2b=fft(vVs2b);
% Vs2b=Vs2b(1:(samples/2+1),1);
% Vs2b=abs(Vs2b)/(samples/2);
% 
% Vs2c=fft(vVs2c);
% Vs2c=Vs2c(1:(samples/2+1),1);
% Vs2c=abs(Vs2c)/(samples/2);
% 
% Vr1a=fft(vVr1a);
% Vr1a=Vr1a(1:(samples/2+1),1);
% Vr1a=abs(Vr1a)/(samples/2);
% 
% Vr1b=fft(vVr1b);
% Vr1b=Vr1b(1:(samples/2+1),1);
% Vr1b=abs(Vr1b)/(samples/2);
% 
% Vr1c=fft(vVr1c);
% Vr1c=Vr1c(1:(samples/2+1),1);
% Vr1c=abs(Vr1c)/(samples/2);
% 
% Vr2a=fft(vVr2a);
% Vr2a=Vr2a(1:(samples/2+1),1);
% Vr2a=abs(Vr2a)/(samples/2);
% 
% Vr2b=fft(vVr2b);
% Vr2b=Vr2b(1:(samples/2+1),1);
% Vr2b=abs(Vr2b)/(samples/2);
% 
% Vr2c=fft(vVr2c);
% Vr2c=Vr2c(1:(samples/2+1),1);
% Vr2c=abs(Vr2c)/(samples/2);