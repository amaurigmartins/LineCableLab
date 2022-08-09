
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: ex3_Y.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.

clear all

load ex3_Y    %--> s, SER

%================================================
%=           Passivity Enforcement              =
%================================================ 
opts.Niter_in=2;
opts.parametertype='Y';
opts.plot.s_pass=2*pi*i*linspace(0,3e4,1001).'; 
opts.plot.ylim=[-2e-3 2e-3];
opts.outputlevel=1;
[SER2,bigYfit_passive,opts2]=RPdriver(SER,s,opts);


%=================================================
%= Comparing original model with perturbed model =
%=================================================
s=2*pi*i*linspace(0,3e4,501);Ns=length(s); %New frequency band
for k=1:Ns
  bigYfit(:,:,k)        =fitcalcPRE(s(k),SER.poles,SER.R,SER.D,SER.E);
  bigYfit_passive(:,:,k)=fitcalcPRE(s(k),SER2.poles,SER2.R,SER2.D,SER2.E);
end
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


