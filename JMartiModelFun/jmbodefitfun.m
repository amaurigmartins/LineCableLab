close all
clear all
clc

fname='carson_cp';
load([fname '.mat'])
addpath('bodefit')

line_length = 1000;
f_choice=200000;

% General parameters
Ns = length(f); % Number of samples
% f = logspace(-2,10,Ns); % Frequency (Hz)
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency

[tau0, tau_mps, tau_opt]=findoptimtau(f,g_dis,line_length,ord);

% Zc fit
zconst=Zch_mod(freq_siz,:);
tol = 3; % Tolerance in decibels to set a new pole and/or new zero

for m=1:ord
    fun=abs(Zch_mod(:,m));
    
    % Bode fit
    [P,Z,k] = Bode_process(fun,f,Ns,tol); % Bode subroutine
    as = k.*poly(Z); bs = poly(P); % Polynomials
    [r,p,ks] = residue(as,bs); % Poles, residues and constant term
    TF=isempty(ks);if(TF==1);ks=0;end

    ffit = zeros(1,Ns)';
    for k = 1:length(p)
        ffit = ffit + (r(k)./(s.' - p(k))).';
    end
    ffit = ffit + ks;
    error = abs(fun - ffit);
   
    fitZcOHLT(m).mode = m;
    fitZcOHLT(m).NORD = length(p);
    fitZcOHLT(m).zInf = abs(zconst(m));
    fitZcOHLT(m).pol = -p;
    fitZcOHLT(m).res = r;
    fitZcOHLT(m).ks = ks;
    
    figure(1)
    semilogx(f,abs(fun), 'DisplayName', ['mode #' num2str(m)]);hold all;
    semilogx(f,abs(ffit), 'o', 'DisplayName', ['fit mode #' num2str(m)]);hold all;
end
figure(1)
axis tight
xlabel('Frequency [Hz]')
ylabel('Z_{ch} [\Omega]')
grid on
legend


% A1 fit
A1min=.1;
tol = 3;
tau=tau_opt;
for o=1:ord
    H(:,o)=exp(-g_dis(:,o).*line_length);
end
[ff, H] = fitnextrap(f,H,ord, 'decdecay');
ss=1j.*2.*pi.*ff;

for m=1:ord
    A1 = H(:,m).*exp(1i*2*pi.*ff.*tau(m));
    fun=abs(A1);
    % Bode fit
    [P,Z,k] = Bode_process(fun(fun>=A1min),ff(fun>=A1min),length(ff(fun>=A1min)),tol); % Bode subroutine
    as = k.*poly(Z); bs = poly(P); % Polynomials
    [r,p,ks] = residue(as,bs); % Poles, residues and constant term
    TF=isempty(ks);if(TF==1);ks=0;end
    ffit = zeros(1,length(ff))';
    for k = 1:length(p)
        ffit = ffit + (r(k)./(ss.' - p(k))).';
    end
    ffit = ffit + ks; 
    
    error = abs(fun(1:Ns) - ffit(1:Ns));
    
    fitA1OHLT(m).mode = m;
    fitA1OHLT(m).NORD = length(p);
    fitA1OHLT(m).tauInf = tau_opt(m);
    fitA1OHLT(m).pol = -p;
    fitA1OHLT(m).res = r;
    fitA1OHLT(m).ks = ks;
    
    figure(2)
    semilogx(ff,abs(fun), 'DisplayName', ['mode #' num2str(m)]);hold all;
    semilogx(ff,abs(ffit), 'o', 'DisplayName', ['fit mode #' num2str(m)]);hold all;
end
figure(2)
axis tight
xlabel('Frequency [Hz]')
ylabel('Propagation factor [unitless]')
grid on
legend;

for o=1:ord
    Ti(o,:)=Ti_dis(find(f==f_choice),(o-1)*ord+1:o*ord);
end

pchfname = ['E:\Users\Amauri\Documents\ATPdata\project\Usp\' fname '200kopt.pch'];
% fcontent = punchJMartiCard(ord, fitZcOHLT, fitA1OHLT, Ti, pchfname);

for m=1:ord
    [b,a] = residue(fitZcOHLT(m).res, fitZcOHLT(m).pol, fitZcOHLT(m).zInf);
    s=1i*2*pi*f;
    for k=1:length(s)
        h(k) = sum(fitZcOHLT(m).res ./ (s(k) + fitZcOHLT(m).pol)) + fitZcOHLT(m).zInf;
    end
    
    figure(3)
    subplot(2,1,1)
    semilogx(f,abs(h));hold all
    xlabel('Frequency [Hz]');
    ylabel('Zc magnitude');
    axis tight
    grid on
    
    
    subplot(2,1,2)
    semilogx(f,unwrap(angle(h))*180/pi, 'DisplayName', ['mode #' num2str(m) ' - OHLT']);hold all
    xlabel('Frequency [Hz]');
    ylabel('Zc angle [deg]');
    axis tight
    grid on
    legend;
    
    for k=1:length(s)
        h(k) = sum(fitA1OHLT(m).res ./ (s(k) + fitA1OHLT(m).pol));
    end
    
    figure(4)
    subplot(2,1,1)
    semilogx(f,abs(h));hold all
    xlabel('Frequency [Hz]');
    ylabel('A1 magnitude');
    axis tight
    grid on
    
    
    subplot(2,1,2)
    semilogx(f,unwrap(angle(h))*180/pi,'DisplayName', ['mode #' num2str(m) ' - OHLT']);hold all
    xlabel('Frequency [Hz]');
    ylabel('A1 angle [deg]');
    axis tight
    grid on
    legend;
end