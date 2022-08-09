clear
close all

f = 60;

%% A simple geometry
hi=8;
hj=8;
x=1.5+0.8; %horizontal spacing

%% Ground parameters
rho=100;
mu_rg = 1;
e_rg = 10;

%% Petterson's Formula

Zm_Petterson = pettersonMutualImpedance(f,hi,hj,x,rho,e_rg);
Ym_Petterson = pettersonMutualAdmittance(f,hi,hj,x,rho,e_rg);

ri = 1.25705e-2;

Zs_Petterson = pettersonSelfImpedance(f,hi,ri,rho,e_rg);
Ys_Petterson = pettersonSelfAdmittance(f,hi,ri,rho,e_rg);

%% Comparison

k=1;
for f=1e3:1e2:1e6
    Ym_Wise(k) = wiseMutualAdmittance(f,hi,hj,x,rho,mu_rg,e_rg);
    Ym_Petterson(k) = pettersonMutualAdmittance(f,hi,hj,x,rho,e_rg);
    Zm_Wise(k) = wiseMutualImpedance(f,hi,hj,x,rho,mu_rg,e_rg);
    Zm_Petterson(k) = pettersonMutualImpedance(f,hi,hj,x,rho,e_rg);
    k=k+1;
end

imagErrorZ = (abs(imag(Zm_Petterson) - imag(Zm_Wise))./(abs(imag(Zm_Wise)))).*100;
realErrorZ = (abs(real(Zm_Petterson) - real(Zm_Wise))./(abs(real(Zm_Wise)))).*100;
ErrorY = (abs(Ym_Petterson - Ym_Wise)./(abs(Ym_Wise))).*100;

f=1e3:1e2:1e6;

figure;
semilogx(f,abs(Zm_Wise),f,abs(Zm_Petterson))
ylabel('|Z_{i,j}|')
xlabel('Frequency [Hz]')
legend('Wise','Petterson')
grid on

figure;
semilogx(f,abs(Ym_Wise),f,abs(Ym_Petterson))
ylabel('|Y_{i,j}|')
xlabel('Frequency [Hz]')
legend('Wise','Petterson')
grid on

figure;
semilogx(f,realErrorZ)
ylabel('Diff_{real(Z)}[%]')
xlabel('Frequency [Hz]')
xlim([1e3 1e6])
grid on

figure;
semilogx(f,imagErrorZ)
ylabel('Diff_{imag(Z)}[%]')
xlabel('Frequency [Hz]')
xlim([1e3 1e6])
grid on

% figure;
% semilogx(f,abs(ErrorY))
% ylabel('Diff_{abs(Y)}[%]')
% xlabel('Frequency [Hz]')
% xlim([1e3 1e6])
% grid on
% 
% figure;
% semilogx(f,angle(ErrorY))
% ylabel('Diff_{angle(Y)}[%]')
% xlabel('Frequency [Hz]')
% xlim([1e3 1e6])
% grid on