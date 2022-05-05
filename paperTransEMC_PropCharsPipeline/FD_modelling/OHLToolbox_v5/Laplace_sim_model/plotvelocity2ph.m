format long;

velocity_matlab1=(2*pi*freq)./imag(gamma_dis(:,1));
velocity_matlab2=(2*pi*freq)./imag(gamma_dis(:,2));
velocity_matlab3=(2*pi*freq)./imag(gamma_dis(:,3));
velocity_matlab4=(2*pi*freq)./imag(gamma_dis(:,4));
velocity_matlab5=(2*pi*freq)./imag(gamma_dis(:,5));
velocity_matlab6=(2*pi*freq)./imag(gamma_dis(:,6));
% velocity_matlab7=(2*pi*freq)./imag(gamma(:,7));
% velocity_matlab8=(2*pi*freq)./imag(gamma(:,8));
%figure(1)
%semilogx(freq,velocity_matlab)

%traveltime=length./velocity;
%figure(2)
%semilogx(freq,traveltime)

% eo=8.854187817e-12;
% mo=4*pi*1e-7;
% 
% er_aera=1;
% mr_aera=1;
% s_aera=0;
% 
% er_ghs=1;
% mr_ghs=1;
% s_ghs=0.001;
% 
% % er_mon1=2.71952;
% % mr_mon1=1;
% % s_mon1=0;
% % 
% % er_mon2=6; 
% % mr_mon2=1;
% % s_mon2=0;
% 
% gamma_aera=sqrt(1i*2*pi*freq*mo*mr_aera.*(s_aera+1i*2*pi*freq*eo*er_aera));
% gamma_ghs=sqrt(1i*2*pi*freq*mo*mr_ghs.*(s_ghs+1i*2*pi*freq*eo*er_ghs));
% % gamma_mon1=sqrt(1i*2*pi*freq*mo*mr_mon1.*(s_mon1+1i*2*pi*freq*eo*er_mon1));
% % gamma_mon2=sqrt(1i*2*pi*freq*mo*mr_mon2.*(s_mon2+1i*2*pi*freq*eo*er_mon2));
% 
% velocity_aera=(2*pi*freq)./imag(gamma_aera);
% velocity_ghs=(2*pi*freq)./imag(gamma_ghs);
% velocity_mon1=(2*pi*freq)./imag(gamma_mon1);
% velocity_mon2=(2*pi*freq)./imag(gamma_mon2);

figure(1)
semilogx(freq,velocity_matlab1,freq,velocity_matlab2,freq,velocity_matlab3,freq,velocity_matlab4,freq,velocity_matlab5,freq,velocity_matlab6)