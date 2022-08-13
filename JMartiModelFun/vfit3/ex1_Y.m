% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: ex1_Y.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
clear all

Ns=501;    %Number of frequency samples
Nc=2;      %Size of Y (after reduction)
bigY=zeros(Nc,Nc,Ns);
Y=zeros(4,4);
s=2*pi*i*logspace(1,5,Ns);

%Component values: 
R1=1;    L1=1e-3; C1=1e-6;
R2=5;    L2=5e-3;
R3=1;    C3=1e-6;
L4=1e-3;
R4=1e-2; L5=20e-3;
R6=10;   C6=10e-6;
R7=1;    C7=2e-6;

%Building Y, reduction:
for k=1:Ns
  sk=s(k);
  y1=1/( R1+sk*L1+1/(sk*C1) );
  y2=1/( R2+sk*L2 );
  y3=1/( R3+1/(sk*C3));
  y4=1/( R4+sk*L4);
  y5=1/(sk*L5);
  y6=1/( R6+1/(sk*C6) );
  y7=1/( R7+1/(sk*C7) );
  
  Y(1,1)= y1+y3;
  Y(2,2)= y4;
  Y(3,3)= y3 +y4 +y5 +y6;
  Y(4,4)= y1 +y2 +y6 +y7;
  
  Y(1,3)=-y3; Y(1,4)=-y1;
  Y(2,3)=-y4;
  Y(3,1)=-y3; Y(3,2)=-y4; Y(3,4)=-y6;
  Y(4,1)=-y1; Y(4,3)=-y6;
  
  %Eliminating nodes 3 and 4:
  Yred=Y(1:2,1:2)-Y(1:2,3:4)*Y(3:4,3:4)^(-1)*Y(3:4,1:2);
  bigY(:,:,k)=Yred;
  %bigY(:,:,k)=diag(diag(Yred));
end  

%================================================
%=           POLE-RESIDUE FITTING               =
%================================================ 
opts.N=8;              %Order of approximation. (Is used when opts.poles=[]).  
opts.poletype='logcmplx';   %Will use logarithmically spaced, complex poles. (Is used when opts.poles=[]).
poles=[]; %[] -->N initial poles are automatically generated as defined by opts.startpoleflag 
[SER,rmserr,bigYfit,opts2]=VFdriver(bigY,s,poles,opts); %Creating state-space model and pole-residue model 

%================================================
%=           Generating external model          =
%================================================ 
NOD='A'; 
fname ='RLC_ATP.txt';
netgen_ATP(SER,NOD,fname);  %Creating branch-cards for ATP 

% Plotting step response: 
NN=length(SER.A) ; I=ones(NN,1);
t=(0:1e-5:5e-3); Nt=length(t);
for k=1:Nt
  if opts2.cmplx_ss==1  
    y=SER.C*diag( diag(SER.A).^(-1).*(exp(diag(SER.A).*t(k))-I) )*SER.B +SER.D;
  else  
    y=SER.C*( (SER.A)^(-1))*((expm((full(SER.A)).*t(k))-diag(I)) ) *SER.B +SER.D;  
  end  
  yy(k)=y(2,1);
end  
figure(4),plot(1000*t,yy);      
xlabel('Time [ms]'); ylabel('Current [A]'); 


