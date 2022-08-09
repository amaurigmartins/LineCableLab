format long;

velocity_matlab=(2*pi*freq)./imag(gamma);
%figure(1)
%semilogx(freq,velocity_matlab)

%traveltime=length./velocity;
%figure(2)
%semilogx(freq,traveltime)

eo=8.854e-12;
mo=12.566370614359172e-7;

er_aera=1;
mr_aera=1;
s_aera=0;

er_ghs=10;
mr_ghs=1;
s_ghs=0.001;

er_mon1=8;
mr_mon1=1;
s_mon1=0;

gamma_aera=sqrt(1i*2*pi*freq*mo*mr_aera.*(s_aera+1i*2*pi*freq*eo*er_aera));
gamma_ghs=sqrt(1i*2*pi*freq*mo*mr_ghs.*(s_ghs+1i*2*pi*freq*eo*er_ghs));
gamma_mon1=sqrt(1i*2*pi*freq*mo*mr_mon1.*(s_mon1+1i*2*pi*freq*eo*er_mon1));

velocity_aera=(2*pi*freq)./imag(gamma_aera);
velocity_ghs=(2*pi*freq)./imag(gamma_ghs);
velocity_mon1=(2*pi*freq)./imag(gamma_mon1);

figure(1)
semilogx(freq,velocity_matlab,freq,velocity_aera,freq,velocity_ghs,freq,velocity_mon1)