function [Y]=calc_Y_Ametani_3ph(freq)

if freq==0
    freq=0.000001;
end

e0=8.854e-12;
er_ins1=2.71952;
er_ins2=6;

r2=0.0119;
r3=0.01855;
r4=0.01975;
r5=0.02275;
h=r5;

D12=2*r5;
D12_icon=sqrt((D12^2)+(2*h)^2);

D13=2*D12;
D13_icon=sqrt((D13^2)+(2*h)^2);

D23=2*r5;
D23_icon=sqrt((D23^2)+(2*h)^2);

Pc=(1/(2*pi*e0*er_ins1))*log(r3/r2);
Ps=(1/(2*pi*e0*er_ins2))*log(r5/r4);

Pi1=[Pc+Ps Ps;Ps Ps;];
Pi2=Pi1;
Pi3=Pi1;

Pi=[Pi1 zeros(2,2) zeros(2,2);zeros(2,2) Pi2 zeros(2,2);zeros(2,2) zeros(2,2) Pi3];

Po_11=(1/(2*pi*e0))*log((2*h)/r5);
Po_12=(1/(2*pi*e0))*log(D12_icon/D12);
Po_13=(1/(2*pi*e0))*log(D13_icon/D13);
Po_23=(1/(2*pi*e0))*log(D23_icon/D23);

Po11=[Po_11 Po_11;Po_11 Po_11;];
Po22=Po11;
Po33=Po11;

Po12=[Po_12 Po_12;Po_12 Po_12;];
Po13=[Po_13 Po_13;Po_13 Po_13;];
Po23=[Po_23 Po_23;Po_23 Po_23;];

Po=[Po11 Po12 Po13;Po12 Po22 Po23;Po13 Po23 Po33];

P=Pi+Po;

Y=1j*2*pi*freq*inv(P);