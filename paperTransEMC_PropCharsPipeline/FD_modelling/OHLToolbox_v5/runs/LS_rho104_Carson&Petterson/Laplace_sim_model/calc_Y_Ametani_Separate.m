e0=8.854e-12;
er_ins1=3.5;
er_ins2=8;

r2=0.0234;
r3=0.0385;
r4=0.0413;
r5=0.0484;
h=r5;

Pc=(1/(2*pi*e0*er_ins1))*log(r3/r2);
Ps=(1/(2*pi*e0*er_ins2))*log(r5/r4);
Po=(1/(2*pi*e0))*log((2*h)/r5);

P=[Pc+Ps+Po Ps+Po; Ps+Po Ps+Po];

num_files=587;
ord=2;
freq=importdata('freq.txt');

Y_tot=zeros(num_files,ord^2);
for o=1:1:num_files
    Y=1j*2*pi*freq(o)*inv(P);
    
    for q=1:1:ord
    Y_dis(1,(q-1)*ord+1:q*ord)=Y(q,:);
    end
    
    Y_tot(o,:)=Y_dis;
end