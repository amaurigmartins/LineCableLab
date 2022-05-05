clear
close all

%% A simple geometry
hi=20;
hj=20;
x=50; %horizontal spacing

%% Ground parameters
rho=1000;
mu_rg = 1;
e_rg = 10;
e0 = 8.854187817e-12;


%% Frequency range 
points_per_dec=10;
f_dec1=1:10/points_per_dec:10;
f_dec2=10:100/points_per_dec:100;
f_dec3=100:1000/points_per_dec:1000;
f_dec4=1000:10000/points_per_dec:10000;
f_dec5=10000:100000/points_per_dec:100000;
f_dec6=100000:1000000/points_per_dec:1000000;
f_range=([f_dec1(1:length(f_dec1)-1) f_dec2(1:length(f_dec2)-1) f_dec3(1:length(f_dec3)-1) f_dec4(1:length(f_dec4)-1) f_dec5(1:length(f_dec5)-1) f_dec6]);

%% Compute stuff
h=[hi hj];
d=[0 x; x 0];
r=[0.01257; 0.01257];
ord=2;

for k=1:1:length(f_range)
    f=f_range(k);
    omega=2*pi*f;
    
    Ym_Wise(k) = wiseMutualAdmittance(f,hi,hj,x,rho,mu_rg,e_rg);
    Zm_Wise(k) = wiseMutualImpedance(f,hi,hj,x,rho,mu_rg,e_rg);
    Zm_Carson(k) = carsonmutual(hi,hj,x,f,1/rho,1);
    Ym_Image(k) = imageMutualAdmittance(f,hi,hj,x);
    
    Zm=Z_pet_mut(h,d,e_rg,mu_rg,1/rho,omega,ord);
    Zpg=Z_self_mut_pg(h,d,r,omega,ord);  
    Z_temp=Zm+Zpg;
    Zm_peterss(k)=Z_temp(1,2);
    
    
    Pm_pet_perf=P_pet_mut_perf(h,d,2);
    Pm_pet_imperf=P_pet_mut_imperf(h,d,e_rg,mu_rg,1/rho,omega,ord);   
    Pm_pet=(1/(2*pi*e0))*(Pm_pet_perf+Pm_pet_imperf);
    Ym_temp = 1i*omega*(inv(Pm_pet));
    Ym_peterss(k)=Ym_temp(1,2);

end

load('Zm_pet_test_50.mat');
load('Ym_pet_test_50.mat');

figure;
semilogx(f_range,abs(Zm_Wise),f_range,abs(Zm_Carson),f_range,abs(Zm_peterss),f_range,abs(Zm_pet))
ylabel('|Z_{i,j}|')
xlabel('Frequency [Hz]')
legend('Wise (NumInt)','Carson','Petersson', 'OHTL Ztot')
grid on

figure;
semilogx(f_range,abs(Ym_Wise),f_range,abs(Ym_Image),f_range,abs(Ym_peterss),f_range,abs(Ym_pet))
ylabel('|Y_{i,j}|')
xlabel('Frequency [Hz]')
legend('Wise (NumInt)','Image','Petersson', 'OHTL Ytot')
grid on

