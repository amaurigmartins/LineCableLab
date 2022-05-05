function [Y]=calc_Y_Ametani(freq)

e0=8.854e-12;
er_ins1=2.71952;
er_ins2=6;

r2=0.0119;
r3=0.01855;
r4=0.01975;
r5=0.02275;
h=r5;

Pc=(1/(2*pi*e0*er_ins1))*log(r3/r2);
Ps=(1/(2*pi*e0*er_ins2))*log(r5/r4);
Po=(1/(2*pi*e0))*log((2*h)/r5);

P=[Pc+Ps+Po Ps+Po; Ps+Po Ps+Po];
Y=1j*2*pi*freq*inv(P);

% Y_tot=zeros(num_files,ord^2);
% for o=1:1:num_files
%     Y=1j*2*pi*freq(o)*inv(P);
%     
%     for q=1:1:ord
%     Y_dis(1,(q-1)*ord+1:q*ord)=Y(q,:);
%     end
%     
%     Y_tot(o,:)=Y_dis;
% end