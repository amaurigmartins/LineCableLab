function [vel_aera,vel_ghs,vel_mon1,vel_mon2]=calc_media_SC_cable(freq)

eo=8.854187817e-12;
mo=4*pi*1e-7;

er_aera=1;
mr_aera=1;
s_aera=0;

er_ghs=1;
mr_ghs=1;
s_ghs=0.01;

er_mon1=2.71952;
%er_mon1=2.3;
mr_mon1=1;
s_mon1=0;

er_mon2=6; 
mr_mon2=1;
s_mon2=0;

gamma_aera=sqrt(1i*2*pi*freq*mo*mr_aera.*(s_aera+1i*2*pi*freq*eo*er_aera));
gamma_ghs=sqrt(1i*2*pi*freq*mo*mr_ghs.*(s_ghs+1i*2*pi*freq*eo*er_ghs));
gamma_mon1=sqrt(1i*2*pi*freq*mo*mr_mon1.*(s_mon1+1i*2*pi*freq*eo*er_mon1));
gamma_mon2=sqrt(1i*2*pi*freq*mo*mr_mon2.*(s_mon2+1i*2*pi*freq*eo*er_mon2));

vel_aera=(2*pi*freq)./imag(gamma_aera);
vel_ghs=(2*pi*freq)./imag(gamma_ghs);
vel_mon1=(2*pi*freq)./imag(gamma_mon1);
vel_mon2=(2*pi*freq)./imag(gamma_mon2);