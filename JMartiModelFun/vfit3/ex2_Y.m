
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: ex2_Y.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.


clear all
load ex2_Y %-->s, bigY


%================================================
%=           POLE-RESIDUE FITTING               =
%================================================ 
opts.N=50 ;%           %Order of approximation. 
opts.poletype='linlogcmplx'; %Mix of linearly spaced and logarithmically spaced poles
opts.weightparam=5; %5 --> weighting with inverse magnitude norm
opts.Niter1=7;    %Number of iterations for fitting sum of elements (fast!) 
opts.Niter2=4;    %Number of iterations for matrix fitting 
opts.asymp=2;      %Fitting includes D   
opts.logx=0;       %=0 --> Plotting is done using linear abscissa axis 
poles=[];      
[SER,rmserr,bigYfit,opts2]=VFdriver(bigY,s,poles,opts); %Creating state-space model and pole-residue model 


%================================================
%=           Passivity Enforcement              =
%================================================ 
clear opts;
opts.parametertype='Y';
opts.plot.s_pass=2*pi*i*linspace(0,2e5,1001).'; 
opts.plot.ylim=[-2e-3 2e-3];

[SER,bigYfit_passive,opts3]=RPdriver(SER,s,opts);


%=================================================
%= Comparing original model with perturbed model =
%=================================================
figure(11),
Nc=length(SER.D);
for row=1:Nc
  for col=row:Nc  
    dum1=squeeze(bigYfit(row,col,:));
    dum2=squeeze(bigYfit_passive(row,col,:));
    h1=semilogy(s/(2*pi*i),abs(dum1),'b'); hold on
    h2=semilogy(s/(2*pi*i),abs(dum2),'r--');     
    h3=semilogy(s/(2*pi*i),abs(dum2-dum1),'g-');         
  end
end 
hold off
xlabel('Frequency [Hz]'); ylabel('Admittance [S]');
legend([h1 h2 h3],'Original model','Perturbed model','Deviation');
