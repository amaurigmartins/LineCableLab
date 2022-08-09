%================================================
% Main program of BODE-DGN
% Authors: Eduardo Salvador Bañuelos Cabral
% José Alberto Gutiérrez Robles
% Bjørn Gustavsen
%===============================================
clc
clear all
close all
%% Initial Settings
Ns = 400; % Number of samples
f = logspace(-2,10,Ns); % Frequency (Hz)
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency
%% Synthetic functions
choice = menu('CHOOSE A FUNCTION','Function F1(s)', 'Function F2(s)','Function F3(s)','Function F4(s)','Function F5(s)','Function F6(s)');
if choice == 1
Fs = (110.*s./((s+10).*(s+100)));
tol = 2.2; % Tolerance in decibels to set a new pole and/or new zero
ite = 40; % Number of iteration in DGN method
end
if choice == 2
Fs = ((s+10).*(s+100))./((s+14580).*(s+10355550));
tol = 2.2; % Tolerance in decibels to set a new pole and/or new zero
ite = 40; % Number of iteration in DGN method
end
if choice == 3
Fs = ((s+10655).*(s+10))./((s+148450).*(s+198545852).*(s+155222220));
tol = 2.2; % Tolerance in decibels to set a new pole and/or new zero
ite = 40; % Number of iteration in DGN method
end
if choice == 4
Fs = (10124).*s.*(s+35.7).*(s+88.9)./((s+100.5).*(s+220.7).*(s+5900).*(s+1370).*(s+21000));
tol = 2; % Tolerance in decibels to set a new pole and/or new zero
ite = 50; % Number of iteration in DGN method
end
if choice == 5
Fs = ((s+79).*(s+1045))./((s+1458).*(s+103555).*(s+127710355).*(s+1244103555));
tol = 1.8; % Tolerance in decibels to set a new pole and/or new zero
ite = 40; % Number of iteration in DGN method
end
if choice == 6
Fs = ((s+64518).*(s+8451629).*(s+312).*(s+54841216192))./((s+456).*(s+7852).* (s+982365).* (s+93256888).*(s+79325684588536));
tol = 2; % Tolerance in decibels to set a new pole and/or new zero
ite = 40; % Number of iteration in DGN method
end
%% Bode method
[P,Z,k] = Bode_process(Fs,f,Ns,tol); % Bode subroutine
as = k.*poly(Z); bs = poly(P); % Polynomials
[r,p,ks] = residue(as,bs); % Poles, residues and constant term
TF=isempty(ks);if(TF==1);ks=0;end
%% Damped Gauss-Newton
Xmcp = [r; p; ks]; % Poles, residues and constant term from Bode
Np = length(p); % Number of poles
[KGN,RGN,PGN,Ff] = Damped_Gauss_Newton(Np,Ns,Fs,s,Xmcp,ite);
Fs_fitDGN = zeros(1,Ns);
for k = 1:length(PGN)
Fs_fitDGN = Fs_fitDGN + (RGN(k)./(s.' - PGN(k))).';
end
Fs_fitDGN = Fs_fitDGN + KGN;
error = abs(Fs - Fs_fitDGN);
%% Plots
figure(3)
loglog(f,abs(Fs),'k',f,abs(Fs_fitDGN),'r--',f,error,'--b','LineWidth',2);
xlabel('Frequency [Hz]'); ylabel('Magnitude [p.u.]');
legend('Data','DGN','Deviation')
figure (4)
semilogy(Ff,'--b','LineWidth',2)
xlabel('Iteration count'); ylabel('RMS-error');