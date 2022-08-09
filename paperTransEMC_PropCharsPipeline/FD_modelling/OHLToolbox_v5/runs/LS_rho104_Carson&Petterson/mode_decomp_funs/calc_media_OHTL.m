function [vel_aera,vel_ghs]=calc_media_OHTL(freq)

eo=8.854187817e-12;
mo=4*pi*1e-7;

er_aera=1;
mr_aera=1;
s_aera=0;

er_ghs=10;
mr_ghs=1;
s_ghs=0.01;

gamma_aera=sqrt(1i*2*pi*freq*mo*mr_aera.*(s_aera+1i*2*pi*freq*eo*er_aera));
gamma_ghs=sqrt(1i*2*pi*freq*mo*mr_ghs.*(s_ghs+1i*2*pi*freq*eo*er_ghs));

vel_aera=(2*pi*freq)./imag(gamma_aera);
vel_ghs=(2*pi*freq)./imag(gamma_ghs);