close all
clear
clc

fname='carson_cp';
load([fname '.mat'])

addpath('vfit3')

line_length = 1000;
f_choice=200000;
Ns=length(f);

% General flags for VF
opts.stable=0;
opts.asymp=2; %1-> d=e=0, 2 -> d!=0, e=0, 3-> d!=0, e!=0
opts.spy2=0;
opts.replace_complex_pairs=0; % Replace complex poles as in https://doi.org/10.1007/s00202-019-00807-8
opts.weightscheme=1;
opts.firstguesstype=1; %real valued poles
opts.output_messages=true;
opts.passive=0;
ERR=1/100;
TOL=2;

% Zc fit
fundb=@(x) 20*log10(abs(x));
zinf=Zch_mod(freq_siz,:);
for m=1:ord
    fun=(Zch_mod(:,m));

    tol=-5;
    NORD=10;
    fit=rationalfit(f,fun,'Tolerance',tol,'TendsToZero',false,'DelayFactor',0,'Npoles',NORD);
    
%     figure(1);semilogx(f,abs(fun),'ko');hold all
    ispassive(fit)
    pol=fit.A;
    p0=remove_complex_poles(pol);
%     fit.A=p0;
    ks=fit.D;
    res=fit.C;
%     [b,a]=residue(fit.C,fit.A,fit.D);
    
%     figure(1);semilogx(f,abs(freqresp(fit,f)));hold all;
    [pol, res, infval, NORD, ffit, err] = vectfit_wrapper(fun,f,ERR,opts);

% [pol, res, infval, NORD, ffit, err]=vf_wrapper(fun,f,p0,1,opts);
%      figure(1);hold all;semilogx(f,abs(ffit))
%     [b,a] = residue(res, pol, infval);
    fitZcOHLT(m).mode = m;
    fitZcOHLT(m).NORD = NORD;
    fitZcOHLT(m).ks = ks;
    fitZcOHLT(m).pol = -p0;
    fitZcOHLT(m).res = res;
    
    figure(1)
    semilogx(f,abs(Zch_mod(:,m)),'o', 'DisplayName', ['mode #' num2str(m)]);hold all;
%     semilogx(f,abs(ffit), 'DisplayName', ['fit mode #' num2str(m)]);hold all;
    semilogx(f,abs(freqresp(fit,f)), 'DisplayName', ['fit mode #' num2str(m)]);hold all;
end
figure(1)
axis tight
xlabel('Frequency [Hz]')
ylabel('Z_{ch} [\Omega]')
grid on
legend


% A1 fit
vel(1:ord)=(2*pi*f(freq_siz))./imag(g_dis(freq_siz,:));
tau(1:ord)=line_length./vel;
asymp_FLAG=1; %d=e=0;
for m=1:ord
    A1 = exp(-g_dis(:,m).*line_length).*exp(1i*2*pi.*f.*tau(m));
    fun=A1;
    [pol, res, infval, NORD, ffit, err]=vectfit_wrapper(fun,f,ERR,opts);
%     [b,a] = residue(res, pol, 0);
    fitA1OHLT(m).mode = m;
    fitA1OHLT(m).NORD = NORD;
    fitA1OHLT(m).tauInf = tau(m);
    fitA1OHLT(m).pol = -pol;
    fitA1OHLT(m).res = res;
    
    figure(2)
    semilogx(f,abs(A1), 'DisplayName', ['mode #' num2str(m)]);hold all;
    semilogx(f,abs(ffit), 'o', 'DisplayName', ['fit mode #' num2str(m)]);hold all;
end
figure(2)
axis tight
xlabel('Frequency [Hz]')
ylabel('Propagation factor [unitless]')
grid on
legend;

for o=1:ord
    Ti(o,:)=Ti_dis(freq_siz,(o-1)*ord+1:o*ord);
end


%%% Punch the PCH card with data from OHTL
% fname = 'E:\Users\Amauri\Documents\ATPdata\projects\Usp\ohtljm.pch';
% fcontent = punchJMartiCard(ord, fitZcOHLT, fitA1OHLT, Ti, fname);


ff=1:1000:100e6;
%%% Plot frequency response from poles and residues to check
for m=1:ord
    % [b,a] = residue(fitZcOHLT(m).res, fitZcOHLT(m).pol, fitZcOHLT(m).zInf);
    s=1i*2*pi*ff;
%     for k=1:length(s)
%         h(k) = sum(fitZcOHLT(m).res ./ (s(k) + fitZcOHLT(m).pol)) + fitZcOHLT(m).zInf;
%     end
%     
%     figure(3)
%     subplot(2,1,1)
%     semilogx(f,abs(h));hold all
%     xlabel('Frequency [Hz]');
%     ylabel('Zc magnitude');
%     axis tight
%     grid on
%     
%     
%     subplot(2,1,2)
%     semilogx(ff,angle(h)*180/pi, 'DisplayName', ['mode #' num2str(m) ' - OHLT']);hold all
%     xlabel('Frequency [Hz]');
%     ylabel('Zc angle [deg]');
%     axis tight
%     grid on
%     legend;
    
%     for k=1:length(s)
%         h(k) = sum(fitA1OHLT(m).res ./ (s(k) + fitA1OHLT(m).pol));
%     end
    
    figure(4)
    subplot(2,1,1)
    semilogx(ff,abs(h));hold all
    xlabel('Frequency [Hz]');
    ylabel('A1 magnitude');
    axis tight
    grid on
    
    
    subplot(2,1,2)
    semilogx(ff,angle(h)*180/pi,'DisplayName', ['mode #' num2str(m) ' - OHLT']);hold all
    xlabel('Frequency [Hz]');
    ylabel('A1 angle [deg]');
    axis tight
    grid on
    legend;
end


%%% Now let's check with data coming from the original ATP file
% rewrite_PCH_data

function [pol, res, ks, N, fit, rmserr] = vf_wrapper(f,freq,poles,N,opts)

addpath('vfit3')

warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:rankDeficientMatrix')

if size(f,1) > 1
    f=f';
end

if size(freq,1) > 1
    freq=freq';
end

Ns=length(freq);
s=1i*2*pi*freq;

weight=ones(1,Ns);

opts.relax=1;
opts.stable=1;
% opts.asymp=1;
opts.spy1=0;
opts.spy2=0;
opts.cmplx_ss=1;

for i=1:N
    [SER,poles,rmserr,fit]=vectfit4(f,s,poles,weight,opts);
end

[R,a]=ss2pr(SER.A,SER.B,SER.C);
res=squeeze(R);
pol=a;
ks = SER.D;

N = length(pol);

warning('on', 'MATLAB:nearlySingularMatrix')
warning('on', 'MATLAB:rankDeficientMatrix')

end