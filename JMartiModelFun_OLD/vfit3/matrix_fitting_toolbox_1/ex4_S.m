
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: ex4_S.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.

clear all
load ex4_S %-->s, bigS


%================================================
%=           POLE-RESIDUE FITTING               =
%================================================ 
opts.N=50 ;%           %Order of approximation. 
opts.poletype='linlogcmplx'; %Mix of linearly spaced and logarithmically spaced poles
opts.Niter1=7;    %Number of iterations for fitting sum of elements (fast!) --> Improved initial poles
opts.Niter2=4;    %Number of iterations for matrix fitting 
opts.asymp=2;      %Fitting includes D   
opts.logx=0;       %=0 --> Plotting is done using linear abscissa axis 
poles=[];      
[SER,rmserr,bigSfit,opts2]=VFdriver(bigS,s,poles,opts); %Creating state-space model and pole-residue model 


%================================================
%=           Passivity Enforcement              =
%================================================ 
clear opts;
opts.plot.s_pass=2*pi*i*linspace(0,2e5,1001).'; 
opts.plot.ylim=[0.95 1.05]; opts.Niter_out=20; 
opts.Niter_in=0;
opts.parametertype='S';
opts.outputlevel=0; %Min. output to screen
[SER,bigSfit_passive,opts3]=RPdriver(SER,s,opts);


%=================================================================
%= Comparing original model with perturbed model in fitting band =
%=================================================================
figure(11),
Nc=length(SER.D);
for row=1:Nc
  for col=row:Nc  
    dum1=squeeze(bigSfit(row,col,:));
    dum2=squeeze(bigSfit_passive(row,col,:));
    h1=semilogy(s/(2*pi*i),abs(dum1),'b'); hold on
    h2=semilogy(s/(2*pi*i),abs(dum2),'r--');     
    h3=semilogy(s/(2*pi*i),abs(dum2-dum1),'g-');         
  end
end 
hold off
xlabel('Frequency [Hz]'); ylabel('Scattering matrix');
legend([h1 h2 h3],'Original model','Perturbed model','Deviation');
