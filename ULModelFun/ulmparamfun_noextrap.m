close all
clear all
clc

addpath('vfit3')

load('matlab.mat','Ych','Ti_dis','g_dis','ord','freq_siz','f','Zmod_src','Ymod_src');
load('matlab.mat',Zmod_src,Ymod_src);

% load('fitulm.mat')
line_length = 1000;

% General parameters
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency


% First build the basic matrices
[H_tmp,P_tmp,H_dis,P_dis,tau]=calc_prop_function(Ti_dis,g_dis,line_length,ord,freq_siz,f);

% for o=1:ord
%    figure(1);semilogx(f,abs(P_tmp(:,o)));hold all
% end

for k=1:freq_siz
    for o=1:ord
        Yc_ph(o,:,k)=Ych(k,(o-1)*ord+1:o*ord);
        Ti(o,:,k)=Ti_dis(k,(o-1)*ord+1:o*ord);
        H_ph(o,:,k)=H_dis(k,(o-1)*ord+1:o*ord);
        H_m(:,o)=exp(-g_dis(:,o).*line_length);
        P_m(:,o)=exp(-g_dis(:,o).*line_length).*exp(1i*2*pi.*f.*tau(o));
    end
    invTi(:,:,k)=inv(Ti(:,:,k));
    trYc(k,1)=trace(Yc_ph(:,:,k));
    H_mod(:,:,k)=diag(H_tmp(k,:));
    P_mod(:,:,k)=diag(P_tmp(k,:));
    trPm(k,1)=trace(P_mod(:,:,k));
end

% Idempotent matrix Dj, eqn. (13) - zanon2021.pdf
for o=1:ord
   for k=1:freq_siz
       mode(o).D(:,:,k)=Ti(:,o,k)*invTi(o,:,k);
       mode(o).H_m(:,:,k)=diag(H_m(k,:));
       mode(o).P_m(:,:,k)=diag(P_m(k,:));
   end
end

% Rewrite  H from idempotent matrix D
H=zeros(ord,ord,freq_siz);
for o=1:ord
    for k=1:freq_siz
        H(:,:,k)=H(:,:,k)+mode(o).D(:,:,k)*diag(exp(-g_dis(k,o).*line_length));
    end
end
test=H-H_ph;
test(abs(test)<0.0001)=0;
any(any(any(test)))

% Find the poles of trYc
tol=-40;
Np=20;
fit=rationalfit(f,trYc,'Tolerance',tol,'TendsToZero',false,'DelayFactor',0,'Npoles',Np);
figure;semilogx(f,abs(freqresp(fit,f)));hold all;semilogx(f,abs(trYc),'ko')
ispassive(fit)
pol=fit.A;

% Now fit all elements of matrix Yc using poles above
opts.N=length(pol);
opts.plot=0;
opts.screen=0;
opts.stable=1;
opts.asymp=2; %D~=0, E=0  ('Proper')
opts.Niter1=1;
opts.Niter2=1;
[SER,rmserr,Yfit,opts2]=VFdriver(Yc_ph,s,pol,opts);
RPopts.Niter_out = 10;
RPopts.Niter_in = 10;
RPopts.parametertype='Y';
RPopts.cmplx_ss=1;
RPopts.outputlevel = 0;
% [SER,~,~]=RPdriver(SER,s,RPopts);

% Write ULM variables for Yc
polesTr=pol;
rYc=SER.R;
resid0Yc=SER.D;

% Find the poles of matrix P_mod (diagonal) - eqn. (2)
% Accurate_transmission_line_modeling_through_optima.pdf
opts.Niter1=100;
opts.Niter2=100;
opts.poletype='logcmplx';
opts.plot=0;
opts.screen=0;
opts.asymp=1; % D=0,  E=0  ('Strictly proper') 
[SER,rmserr,Pfit,opts2]=VFdriver(P_mod,s,[],opts);
RPopts.Niter_out = 10;
RPopts.Niter_in = 10;
RPopts.parametertype='Y';
RPopts.cmplx_ss=1;
RPopts.outputlevel = 1;
[SER,~,~]=RPdriver(SER,s,RPopts);
pol=SER.poles;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % this is to test if fit agrees with original values
for k=1:freq_siz
   trPfit(k,1)=trace(Pfit(:,:,k));
   for o=1:ord
       tmpPfit(k,o)=Pfit(o,o,k);
   end
end
figure(30);semilogx(f,abs(trPfit));hold all; semilogx(f,abs(trPm),'ko')

for o=1:ord
    figure(40);semilogx(f,abs(tmpPfit(:,o)));hold all; semilogx(f,abs(P_tmp(:,o)),'ko')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now find the residues of phase-domain matrix H using poles above 
opts.Niter1=1;
opts.Niter2=1;
for o=1:ord
    clear tmpfun
    for k=1:freq_siz
        tmpfun(:,:,k)=mode(o).D(:,:,k)*P_mod(:,:,k);
    end
    [SER,rmserr,Hphfit,opts2]=VFdriver(tmpfun,s,pol,opts);
    Cij(:,:,:,o)=SER.R;
end

% Write ULM variables for H
pPj=pol;
rA=Cij;
tauOtim=tau;
NFase=ord;

% Punch ULM data file
GeraTXT
