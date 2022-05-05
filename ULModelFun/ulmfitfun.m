close all
clear all
clc

addpath('vfit3')

load('matlab.mat')
% load('fitulm.mat')
line_length = 1000;

% General parameters
Ns = length(f); % Number of samples
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency


% First build the basic matrices
[H_mod,P_mod,H_dis,P_dis,tau]=calc_prop_function(Ti_dis,g_dis,line_length,ord,freq_siz,f);

for k=1:freq_siz
    for o=1:ord
        Yc_ph(o,:,k)=Ych(k,(o-1)*ord+1:o*ord);
        Ti(o,:,k)=Ti_dis(k,(o-1)*ord+1:o*ord);
        H_ph(o,:,k)=H_dis(k,(o-1)*ord+1:o*ord);
    end
    invTi(:,:,k)=inv(Ti(:,:,k));
    trYc(k,1)=trace(Yc_ph(:,:,k));
end

% Modal matrix Dj, eqn. (13) - zanon2021.pdf
for o=1:ord
   for k=1:freq_siz
       mode(o).D(:,:,k)=Ti(:,o,k)*invTi(o,:,k);
   end
end

% Rewrite vectors P_mod to ensure all terms -> 0 as freq -> inf
lastdec=log10(f(end));
new_freq=logspace(lastdec,lastdec+4,500);
new_freq_siz=length(f)+length(new_freq(2:end));
ff=[f; new_freq(2:end)'];
ww = 2.*pi.*ff; % Frequency (rad/seg)
ss = 1j*ww; % Complex Frequency
for o=1:ord
    fun=P_mod(:,o);
    exp_model = @(g,x) (exp(-(x.*g)));
    obj_fun = @(params) norm(exp_model(params(1),f)-fun);
    sol = fminsearch(obj_fun, [fun(end),0,1]);
    g_sol = sol(1);
    fun_extrap = exp_model(g_sol, new_freq);
    mode(o).P=[fun; fun_extrap(2:end)'];
%         figure;semilogx(f,abs(fun),'ko');hold all;semilogx([f; new_freq(2:end)'],abs(mode(o).P))
    M=repmat(mode(o).Dj(:,:,freq_siz),[1 1 new_freq_siz-freq_siz-1]);
    mode(o).Dj(:,:,freq_siz+1:new_freq_siz)=repmat(mode(o).Dj(:,:,freq_siz),[1 1 new_freq_siz-freq_siz]);
end
TOL=1e-6;
P_mod_extrap=reshape([mode.P],[new_freq_siz,ord]);
P_mod_extrap(P_mod_extrap<=TOL)=0; %to prevent underflow issues

for k=1:new_freq_siz
    Pmdiag(:,:,k)=diag(P_mod_extrap(k,:));
end

% save fitulm.mat tau Yc Ti invTi Pmdiag P_mod trYc ord f freq_siz

% Trying to fit the entire matrix Yc
% opts.Niter1=4;
% opts.Niter2=4;
% opts.N=5;
% opts.poletype='logcmplx';
% opts.stable=1;
% opts.cmplx_ss=1;
% opts.plot=0;
% opts.screen=0;
% opts.asymp=1;
% [SER,rmserr,Ycfit,opts2]=VFdriver(Yc,s,[],opts);
% pol=SER.poles;

% Using matlab native function
% tol=-50;
% Np=15;
% fit=rationalfit(f,trYc,'Tolerance',tol,'TendsToZero',false,'DelayFactor',0,'Npoles',Np);
% % figure;plot(f,abs(freqresp(fit,f)));
% % ispassive(fit)
% res=fit.C;
% pol=fit.A;
% const=fit.D;

% Using regular vectfit
spyplot_FLAG=0;
msg_FLAG=false;
ERR=1/100;
stable_FLAG=1; %enforce poles to the left side
passive_FLAG=1;
asymp_FLAG=2; %d!=0,e=0
[pol, res, infval, NORD, ffit, err]=vfit3_wrapper(trYc,f,ERR,asymp_FLAG,stable_FLAG,passive_FLAG,spyplot_FLAG,msg_FLAG);


% Now fit all elements of matrix Yc using pol as initial guess
opts.N=length(pol);
opts.plot=0;
opts.screen=0;
opts.stable=stable_FLAG;
opts.asymp=asymp_FLAG;
opts.Niter1=10;
opts.Niter2=10;
[SER,rmserr,Yfit,opts2]=VFdriver(Yc,s,pol,opts);
RPopts.Niter_out = 100;
RPopts.Niter_in = 100;
RPopts.parametertype='Y';
RPopts.cmplx_ss=1;
RPopts.outputlevel = 1;
[SER,~,~]=RPdriver(SER,s,RPopts);
% pol=SER.poles;

% Write ULM variables for Yc
polesTr=pol;
rYc=SER.R;
resid0Yc=SER.D;

% Find the poles of matrix P_mod (diagonal) - eqn. (2)
% Accurate_transmission_line_modeling_through_optima.pdf
opts.Niter1=10;
opts.Niter2=10;
opts.poletype='logcmplx';
opts.plot=0;
opts.screen=0;
opts.asymp=1; % D=0,  E=0  ('Strictly proper') 
[SER,rmserr,Pfit,opts2]=VFdriver(Pmdiag,ss,[],opts);
pol=SER.poles;

% And, hopefully, the residue matrices for each mode
opts.N=length(pol);
opts.plot=0;
opts.screen=0;
opts.asymp=1;
opts.Niter1=1;
opts.Niter2=1;
for o=1:ord
    clear fun
    for k=1:new_freq_siz
          fun(:,:,k)=mode(o).Dj(:,:,k)*P_mod_extrap(k,o);
%         DP(:,:,k)=D(o).values(:,:,k)*P_mod(k,o); %eqn. (15) - zanon2021.pdf
    end
    [SER,rmserr,DPfit,opts2]=VFdriver(fun,ss,pol,opts);
    Cij(:,:,:,o)=SER.R;
end

% Write ULM variables for H
pPj=pol;
rA=Cij;
tauOtim=tau;
NFase=ord;

% Punch ULM data file
GeraTXT



