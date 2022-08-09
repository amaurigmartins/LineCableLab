format long;

%att_1=(real(gamma(:,1))/8.686)*1e3;
%att_2=(real(gamma(:,2))/8.686)*1e3;
att_1=real(gamma(:,1));
att_2=real(gamma(:,2));
att_3=real(gamma(:,3));

eo=8.854e-12;
mo=12.566370614359172e-7;

er_aera=1;
mr_aera=1;
s_aera=0;
 
er_ghs=1;
mr_ghs=1;
s_ghs=0.001;

% er_mon1=2.71952;
% mr_mon1=1;
% s_mon1=0;
%  
% er_mon2=6;
% mr_mon2=1;
% s_mon2=0;
 
gamma_aera=sqrt(1i*2*pi*freq*mo*mr_aera.*(s_aera+1i*2*pi*freq*eo*er_aera));
gamma_ghs=sqrt(1i*2*pi*freq*mo*mr_ghs.*(s_ghs+1i*2*pi*freq*eo*er_ghs));
% gamma_mon1=sqrt(1i*2*pi*freq*mo*mr_mon1.*(s_mon1+1i*2*pi*freq*eo*er_mon1));
% gamma_mon2=sqrt(1i*2*pi*freq*mo*mr_mon2.*(s_mon2+1i*2*pi*freq*eo*er_mon2));
 
att_aera=real(gamma_aera);
att_ghs=real(gamma_ghs);
% att_mon1=(2*pi*freq)./imag(gamma_mon1);
% att_mon2=(2*pi*freq)./imag(gamma_mon2);

figure(1)
semilogx(freq,att_1,freq,att_2,freq,att_3)