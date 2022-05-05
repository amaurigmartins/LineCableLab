function [Y]=calc_Y_Ametani_1ph(freq)

if freq==0
    freq=0.000001;
end

% e0=8.854e-12;
% er_ins1=2.71952;
% er_ins2=6;
% 
% r2=0.0119;
% r3=0.01855;
% r4=0.01975;
% r5=0.02275;
% h=r5;

e0=8.854e-12;
er_ins1=3.7;
er_ins2=6;

r2=0.00065;
r3=0.0029;
r4=0.002956525;
r5=0.00425;
h=r5;

Pc=(1/(2*pi*e0*er_ins1))*log(r3/r2);
Ps=(1/(2*pi*e0*er_ins2))*log(r5/r4);
Po=(1/(2*pi*e0))*log((2*h)/r5);

P=[Pc+Ps+Po Ps+Po; Ps+Po Ps+Po];
Y=1j*2*pi*freq*inv(P);