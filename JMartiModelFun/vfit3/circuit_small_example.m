% Small example demonstrating the intefacing of rational models with circuit
% simulators using discrete convolution.
% 
% - Two-port Y-parameters are computed for a small electrical circuit
% - The Y-parameters are converted into Z-parameters and S-parameters
% - Each parameter set is fitted with a highly accurate rational model 
%   using Vector Fitting.
% - Each rational model is included in a time domain simulation via a
%   Norton equivalent, assuming a fixed time step length. 
% - The obtained simulation results are compared to th result by a 
%   conventional circuit simulator (PSCAD).
% - Results are also shown for a general transfer function, used for
%   calculating the voltage transfer function between two nodes. 
%
%  Download site:
%  http://www.energy.sintef.no/Produkt/VECTFIT/index.asp 
%
%  ===================================================================
%  If using any of these templates in a scientific work, reference should
%  be made to:
%  B. Gustavsen and H.M.J. De Silva, "Inclusion of rational models in an
%  electromagnetic transients program - Y-parameters, Z-parameters,
%  S-parameters, transfer functions", IEEE Trans. Power Delivery, 
%  vol. 28, no. 2, pp. 1164-1174, April 2013.
%  ====================================================================
%
%  17.03.2013. Bjorn Gustavsen, SINTEF Energy Research, Norway.
%  


clear all

% Update this path to where you have placed the Matrix Fitting Toolbox.
addpath d:\user\mtrx_fitter_new2\auxi



%##########################################################################
%=============================
% Simulation: Y-parametre:
%=============================
%##########################################################################

% A=-2*pi*1000;
% C=-A;
% B=1;
% D=0.8;

clear all

Ns=501;    %Number of frequency samples
Nc=2;      %Size of Y (after reduction)
bigY=zeros(Nc,Nc,Ns);
Y=zeros(4,4);
s=2*pi*i*logspace(1,5,Ns);

%Component values: 
         L1=0.1e-3;
R2=1;    L2=1e-3;   C2=1e-6;
R3=1;               C3=2e-6;
R4=1;               C4=1e-6;
R5=10;              C5=10e-6;
R6=5;    L6=5e-3;
R7=1;               C7=2e-6;
R8=1e-2; L8=1e-3;  
         L9=20e-3;
                    C10=0.2e-6;
  
%Building Y, reduction:
for k=1:Ns
  Y=zeros(5);  
  sk=s(k);
  
  %Branch admittances
  y1 =1/(sk*L1); 
  y2 =1/( R2+sk*L2+1/(sk*C2) );
  y3 =1/( R3+1/(sk*C3) );
  y4 =1/( R4+1/(sk*C4));
  y5 =1/( R5+1/(sk*C5));
  y6 =1/( R6+sk*L6) ;
  y7 =1/( R7+1/(sk*C7) );  
  y8 =1/( R8+sk*L8 );  
  y9 =1/(sk*L9)
  y10=sk*C10;

  %Adding contribution from branch admittances:
  
  Y=add_branch(5,4,y2,Y); 
  Y=add_branch(4,0,y3,Y); 
  Y=add_branch(5,3,y4,Y); 
  Y=add_branch(3,4,y5,Y); 
  Y=add_branch(4,0,y6,Y);  
  Y=add_branch(2,3,y7,Y); 
  Y=add_branch(2,3,y8,Y);  
  Y=add_branch(3,0,y9,Y);  
  Y=add_branch(2,0,y10,Y);   
  Y=add_branch(1,5,y1,Y);  
  
  %Eliminating nodes 3-5:
  Yred=Y(1:2,1:2)-Y(1:2,3:5)*Y(3:5,3:5)^(-1)*Y(3:5,1:2);
  Yred=(Yred+Yred.')/2;
  bigY(:,:,k)=Yred;

end  
  
%================================================
%=           Y: POLE-RESIDUE FITTING               =
%================================================ 
            %Order of approximation. (Is used when opts.poles=[]).  
opts.N=10; %Model order
opts.poletype='logcmplx';  
opts.asymp=3; %Model includes R0 and E
opts.stable=1;
opts.NE=0;
poles=[]; 
[SERY,rmserr,bigYfit,opts2]=VFdriver(bigY,s,poles,opts); 

freq=s/(2*pi*i);
figure(71),
for row=1:Nc
  for col=1:Nc
     dum1=squeeze(bigY(row,col,:)).';
     dum2=squeeze(bigYfit(row,col,:)).';     
     h1=loglog(freq,abs(dum1),'b');hold on
     h2=loglog(freq,abs(dum2),'r--');
     h3=loglog(freq,abs(dum2-dum1),'g-.');
  end
end
hold off
xlabel('Frequency [Hz]')
title('Admittance matrix (magnitude)')
legend([h1(1) h2(1) h3(1)],'Data','Model','Deviation',4)
hold off

%##########################################################################

%================================================
%=           Y: SIMULATION              =
%================================================ 

%+++++++++++++++++++++++++++++++
%**Simulation parameters:
Dt=1e-5; time=(0:Dt:5e-3); Nt=length(time);

%+++++++++++++++++++++++++++++++
%**External circuit:
Rsour=5; 
gsour=1/Rsour;

%+++++++++++++++++++++++++++++++
%**rational model:

%Throw out second part of complex poles:
[SER]=reducecmplx(SERY); 
A=full(diag(SER.A)); %A is now a column vector
B=SER.B
C=SER.C
D=SER.D;
E=SER.E;

%+++++++++++++++++++++++++++++++
%Initialization:
N=length(A);
Nc=length(D);
X=zeros(N,1);

a=A.';
alfa=(1+a*Dt/2)./(1-a*Dt/2);
lamb=(Dt/2)./(1-a*Dt/2);
mu=lamb;
GY=D+real(C*diag(lamb)*B); %Norton conductance for model

dum=alfa.*lamb +mu;
Ctilde=zeros(Nc,N);
for row=1:Nc
  Ctilde(row,:)=C(row,:).*dum; %Ctilde
end  

%Contribution from E-matrix, if present:
if max(max(abs(E)))>0
  alfa=[alfa -ones(1,Nc)];  
  Ctilde=[Ctilde eye(Nc)];
  B=[B;-4*E./Dt];
  GY=GY+2*E./Dt;
  X=[X;zeros(Nc,1)];
end    

GG=zeros(2,2);              %Global conductance matrix
GG(1:2,1:2)=GG(1:2,1:2)+GY; %Contribution from model
GG(1,1)=GG(1,1)+gsour;      %Contribution from external circuit

ZZ=GG^(-1);   %Impedance matrix
Ihis=[0 0].'; %Global history current source


t=0;
for k=1:Nt
  t=t+Dt;

  Iind=[gsour 0].'; %Global vector of independent current sources 1  
  Itot=Iind+Ihis; %Total current
  V=ZZ*Itot;      %Node voltages
  
  %Updating history current source of model
  X=alfa.'.*X +B*V;
  Ihis=-real(Ctilde*X); 
  
  pgbY(1,k)=gsour-gsour*V(1);  %Saving current flowing into node #1
  pgbY(2,k)=V(2);;             %Saving voltage at node #2   

  bigV(k)=V(1); %to be used as excitation in voltage transfer computation
end %for k=1:Nt  


%##########################################################################
%=============================
%  Z-parameters:
%=============================
%##########################################################################

for k=1:Ns
  Y=squeeze(bigY(:,:,k));
  bigZ(:,:,k)=Y^(-1);
end 


%================================================
%=           Z: POLE-RESIDUE FITTING               =
%================================================ 
opts.N=10; %Model order  
poles=[]; 
opts.poletype='logcmplx';   
opts.asymp=3; %Model includes R0 and E
[SERZ,rmserr,bigZfit,opts2]=VFdriver(bigZ,s,poles,opts); 

freq=s/(2*pi*i);
figure(72),
for row=1:Nc
  for col=1:Nc
     dum1=squeeze(bigZ(row,col,:)).';
     dum2=squeeze(bigZfit(row,col,:)).';     
     h1=loglog(freq,abs(dum1),'b');hold on
     h2=loglog(freq,abs(dum2),'r--');
     h3=loglog(freq,abs(dum2-dum1),'g-.');
  end
end
hold off
xlabel('Frequency [Hz]')
title('Impedance matrix (magnitude)')
legend([h1(1) h2(1) h3(1)],'Data','Model','Deviation',4)
hold off


%================================================
%=           Z: SIMULATION              =
%================================================ 
%+++++++++++++++++++++++++++++++
%**Simulation parameters:
Dt=1e-5; time=(0:Dt:5e-3); Nt=length(time);


%+++++++++++++++++++++++++++++++
%**External circuit:
Rsour=5; 
gsour=1/Rsour;

%+++++++++++++++++++++++++++++++
%**rational model:

[SERZ]=reducecmplx(SERZ); 
A=full(diag(SERZ.A));
B=SERZ.B;
C=SERZ.C;
D=SERZ.D;
E=SERZ.E;

%+++++++++++++++++++++++++++++++
%Initialization:
N=length(A);
Nc=length(D);
X=zeros(N,1);

a=A.';
alfa=(1+a*Dt/2)./(1-a*Dt/2);
lamb=(Dt/2)./(1-a*Dt/2);
mu=lamb;

Zthevenin=zeros(Nc);
Zthevenin=real(C*diag(lamb)*B);
Zthevenin=Zthevenin+D;

dum=alfa.*lamb +mu;
Ctilde=zeros(Nc,N);
for row=1:Nc
  Ctilde(row,:)=C(row,:).*dum;
end  
%Contribution from E-matrix, if present:
if max(max(abs(E)))>0
  alfa=[alfa -ones(1,Nc)];  
  Ctilde=[Ctilde eye(Nc)];
  B=[B;-4*E./Dt];
  Zthevenin=Zthevenin+2*E./Dt;
  X=[X;zeros(Nc,1)];
end   
GY=Zthevenin^(-1); %Converting Zthevenin into Gnorton

%Assembling global conductance matrix, calculating impedance matrix:
GG=GY;
GG(1,1)=GG(1,1)+gsour;
ZZ=GG^(-1);   %Impedance matrix
Ihis=[0 0].'; %Global history current source

%Time step loop:
t=0;
V=[0 0].';
for k=1:Nt
  t=t+Dt;
  
  Iind=[gsour 0].'; %Global vector of independent current sources
  Itot=Iind+Ihis;     %Total current
  V=ZZ*Itot;          %Node voltages

  %Updating history currnt source:
  I=(GY*V-Ihis); %current flowing into terminals  
  X=alfa.'.*X +B*I;
  Ihis=GY*real(Ctilde*X); 
  pgbZ(1,k)=I(1);  %saving current flowing into node #1  
  pgbZ(2,k)=V(2);              %Saving voltage at node #2    

end %for k=1:Nt  


%##########################################################################
%=============================
%  S-parameters:
%=============================
%##########################################################################

%================================================
% Konvertering fra Y-parametre til S-parametre:
%================================================

R=diag([100 200]); %Reference impedances
sqR=sqrt(R);       %Square-root of reference impedances
I=eye(Nc); %Nc=2
for k=1:Ns
  Y=bigY(:,:,k);
  bigS(:,:,k)=(I+sqR*Y*sqR)^(-1)*(I-sqR*Y*sqR);
end


%================================================
%=           POLE-RESIDUE FITTING               =
%================================================ 
opts.N=11; %Model order
opts.poletype='logcmplx';  
poles=[]; 
opts.stable=0;
opts.asymp=2; %Model includes R0
[SERS,rmserr,bigSfit,opts2]=VFdriver(bigS,s,poles,opts); 

figure(73),
for row=1:Nc
  for col=1:Nc
     dum1=squeeze(bigS(row,col,:)).';
     dum2=squeeze(bigSfit(row,col,:)).';     
     h1=loglog(freq,abs(dum1),'b');hold on
     h2=loglog(freq,abs(dum2),'r--');
     h3=loglog(freq,abs(dum2-dum1),'g-.');
  end
end
hold off
xlabel('Frequency [Hz]')
title('Scattering matrix (magnitude)')
legend([h1(1) h2(1) h3(1)],'Data','Model','Deviation',4)
hold off


%=============================
% Simulation: S-parametre:
%=============================

%+++++++++++++++++++++++++++++++
%**Rational model:
[SERS]=reducecmplx(SERS); 
A=full(diag((SERS.A)));
B=SERS.B
C=SERS.C
D=SERS.D;

%**Simulation parameters:
Dt=1e-5; time=(0:Dt:5e-3); Nt=length(time);

%+++++++++++++++++++++++++++++++
%**External circuit:
Rsour=5; 
gsour=1/Rsour;

%+++++++++++++++++++++++++++++++
%Initialization:
N=length(A);
Nc=length(D);
X=zeros(N,1);
II=eye(Nc);

a=A.';
alfa=(1+a*Dt/2)./(1-a*Dt/2);
lamb=(Dt/2)./(1-a*Dt/2);
mu=lamb;

G=D+real(C*diag(lamb)*B);
GS= R^(-0.5)*(II-G)*(II+G)^(-1)*R^(-0.5);

dum=alfa.*lamb +mu;
Ctilde=zeros(Nc,N);
for row=1:Nc
  Ctilde(row,:)=C(row,:).*dum;
end  
GAMMA=2*(R^(-0.5))*(II+G)^(-1)*Ctilde;

%+++++++++++++++++++++++++++++++
%Assembling global conductance matrix, calculate impedance matrix:
GG=GS;
GG(1,1)=GG(1,1)+gsour;
ZZ=GG^(-1);
Ihis=[0 0].';

%+++++++++++++++++++++++++++++++
%Time domain simulation:
t=0;
for k=1:Nt
  t=t+Dt;

  Iind=[gsour 0].';
  I=Iind+Ihis;
  V=ZZ*I; %Node voltages
  
  %Calculating incoming power wave
  a=(II+G)^(-1)*( R^(-0.5)*V -real(Ctilde*X) );
  
  %Updating history current source
  X=alfa.'.*X +B*a;
  Ihis=real(GAMMA*X);  
  
  pgbS(1,k)=gsour-gsour*V(1); %Saving currentvflowing into node #1
  pgbS(2,k)=V(2);             %Saving voltage at node #2  

end %for k=1:Nt  


%##########################################################################
%=============================
%  Voltage Transfer function:
%=============================
%##########################################################################

% Voltage transfer function from node 1 to node 2
for k=1:Ns
  Y=squeeze(bigY(:,:,k));  
  bigH(:,:,k)=-(Y(2,2))^(-1)*Y(2,1);
end
%================================================
%=           POLE-RESIDUE FITTING: H            =
%================================================ 
opts.N=11; %Model order  
opts.poletype='logcmplx'; 
poles=[]; 
opts.stable=1
[SERH,rmserr,bigHfit,opts2]=VFdriver(bigH,s,poles,opts); 


figure(74),
for row=1:1
  for col=1:1
     dum1=squeeze(bigH(row,col,:)).';
     dum2=squeeze(bigHfit(row,col,:)).';     
     h1=loglog(freq,abs(dum1),'b');hold on
     h2=loglog(freq,abs(dum2),'r--');
     h3=loglog(freq,abs(dum2-dum1),'g-.');
  end
end
hold off
xlabel('Frequency [Hz]')
title('Voltage transfer function (magnitude)')
legend([h1(1) h2(1) h3(1)],'Data','Model','Deviation',4)
hold off

%================================================
%=           Time domain simulation: H          =
%================================================ 

%+++++++++++++++++++++++++++++++
%**Rational model:
[SERH]=reducecmplx(SERH); 
A=full(diag(SERH.A));
B=SERH.B
C=SERH.C
D=SERH.D;

%**Simulation parameters:
Dt=1e-5; time=(0:Dt:5e-3); Nt=length(time);

%+++++++++++++++++++++++++++++++
%**Initialization:
N=length(A);
Nc=length(D);
X=zeros(N,1);

a=A.';
alfa=(1+a*Dt/2)./(1-a*Dt/2);
lamb=(Dt/2)./(1-a*Dt/2);
mu=lamb;

dum=alfa.*lamb +mu;
Ctilde=zeros(1,N);
for row=1:1
  Ctilde(row,:)=C(row,:).*dum;
end  

GH=real(C*diag(lamb)*B);
GH=GH+D;


%+++++++++++++++++++++++++++++++
%Time domain simulation:
t=0;
oldVmode=zeros(N,1);
for k=1:Nt
  u=bigV(k); %Voltage at node #1 calculated using Y-parameter model 
  t=t+Dt;  
  Vmode=B*u;
  X=alfa.'.*X +oldVmode;
  oldVmode=Vmode;
  y = GH*u +real(Ctilde*X); 
  pgbH(1,k)=0;      
  pgbH(2,k)=y; %The voltage on port #2   
end 



%##########################################################################
%=============================
%  Comparing trime domain results:
%=============================
%##########################################################################

%PSCAD simulation of RLC circuit
load -ASCII C.m;

t_pscad=C(:,1).';  %time 
v_pscad=1e-3*C(:,2).';
i_pscad=1e-3*C(:,3).';

%=============================
% CURRENT at port #1:
%=============================
figure(101),
h1=plot(1e3*time,pgbY(1,:),'b-');hold on
h2=plot(1e3*time,pgbZ(1,:),'r--');
h3=plot(1e3*time,pgbS(1,:),'g-.'); hold off
legend([h1 h2(1) h3(1)],'Y-parameters','Z-parameters','S-parameters',1);
xlabel('Time [ms]');
ylabel('Current [A]')
title('Simulation result') 
hold off


figure(102),
h1=plot(1e3*time,pgbY(1,:)-i_pscad,'b-');hold on
h2=plot(1e3*time,pgbZ(1,:)-i_pscad,'r--');
h3=plot(1e3*time,pgbS(1,:)-i_pscad,'g-.'); hold off
title('Deviation from PSCAD simulation result') 
legend([h1 h2(1) h3(1)],'Y-parameters','Z-parameters','S-parameters',1);
xlabel('Time [ms]');
ylabel('Current [A]')
hold off



%=============================
%  VOLTAGE at port #2:
%=============================
figure(103),
h1=plot(1e3*time,pgbY(2,:),'b-');hold on
h2=plot(1e3*time,pgbZ(2,:),'r--');
h3=plot(1e3*time,pgbS(2,:),'g-.'); 
h4=plot(1e3*time,pgbH(2,:),'k:','linewidth',1.0); hold off
title('Simulation result') 
legend([h1(1) h2(1) h3(1) h4(1)],'Y-parameters','Z-parameters','S-parameters','Transfer function',1);
xlabel('Time [ms]');
ylabel('Voltage [V]')

figure(104),
h1=plot(1e3*time,pgbY(2,:)-v_pscad,'b-');hold on
h2=plot(1e3*time,pgbZ(2,:)-v_pscad,'r--');
h3=plot(1e3*time,pgbS(2,:)-v_pscad,'g-.'); 
h4=plot(1e3*time,pgbH(2,:)-v_pscad,'k:'),hold off
title('Deviation from PSCAD simulation result') 
legend([h1(1) h2(1) h3(1) h4(1)],'Y-parameters','Z-parameters','S-parameters','Transfer function',1);
xlabel('Time [ms]');
ylabel('Voltage [V]')




