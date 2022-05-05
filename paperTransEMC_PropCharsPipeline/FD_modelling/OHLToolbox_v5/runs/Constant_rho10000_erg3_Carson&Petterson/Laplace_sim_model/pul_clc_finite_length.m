% Finite length earth impedance and admittance
function [g,Z,Y,Ti,freq]=pul_clc_finite_length(ord,num_files,x)
%% Geometry data
h1=11.5;
h2=11.5;
h3=11.5;
y12=6.5;
y13=13;
y23=6.5;
r1=0.012573;
r2=0.012573;
r3=0.012573;

%% EM data
m0=4*pi*1e-7;
e0=8.854187817e-12;
res_g=1000;
freq=importdata('freq.txt');
freq(1)=0.000001;

Z=zeros(num_files,ord^2);
Y=zeros(num_files,ord^2);

filedata_Zskin=importdata('fbatchZskin.txt');
path='C:\Program Files\MATLAB\R2010b\work\Multiphase transient model version 3.0 (with nodal method, modal analysis and frequency dependent Ti)';
path_Z=strcat(path,'\Z');

for k=1:1:num_files

    frequency=freq(k);
    omega=2*pi*frequency;
    he=sqrt(res_g./(1i.*omega*m0));

%% Impedance and Admittance Parameters

    [Z11,Po11]=self_impdnce_fin_fun(x,h1,r1,he,omega,m0);
    [Z22,Po22]=self_impdnce_fin_fun(x,h2,r2,he,omega,m0);
    [Z33,Po33]=self_impdnce_fin_fun(x,h3,r3,he,omega,m0);
    [Z12,Po12]=mtl_impdnce_fin_fun(x,h1,h2,y12,he,omega,m0);
    [Z13,Po13]=mtl_impdnce_fin_fun(x,h1,h3,y13,he,omega,m0);
    [Z23,Po23]=mtl_impdnce_fin_fun(x,h2,h3,y23,he,omega,m0);
    
    Zdis=[Z11 Z12 Z13; Z12 Z22 Z23; Z13 Z23 Z33;]/x;
    
    cd(path_Z);
    Ztmp_skin=parse_matrix(k,filedata_Zskin,ord);
    cd(path);
    [Zdis_skin,Ydis]=zyph_param(Ztmp_skin,Ztmp_skin,ord);
    
    Zdis=Zdis+Zdis_skin;
    
    [P11]=self_admttnce_fin_fun(x,h1,r1);
    [P22]=self_admttnce_fin_fun(x,h2,r2);
    [P33]=self_admttnce_fin_fun(x,h3,r3);
    [P12]=mtl_admttnce_fin_fun(x,h1,h2,y12);
    [P13]=mtl_admttnce_fin_fun(x,h1,h3,y13);
    [P23]=mtl_admttnce_fin_fun(x,h2,h3,y23);
    
    Pdis=[P11 P12 P13; P12 P22 P23; P13 P23 P33;];
    Cdis=2*pi*e0.*(x.^2)*inv(Pdis);
    Ydis=1i.*omega*Cdis./x;
          
    for o=1:ord
        Z(k,(o-1)*ord+1:o*ord)=Zdis(o,:);
        Y(k,(o-1)*ord+1:o*ord)=Ydis(o,:);
    end    
end

%[Ti,g]=NR_calc_norm_str(ord,num_files,Z,Y,freq); % Συνάρτηση Newton-Raphson με seeding αυξανομένης της συχνότητας που επιλύει το eigenvector switchover - Υπολογίζει πίνακα ιδιοδιανυσμάτων Ti (current) (num_files x ord^2) και ιδιοτιμές g (sqrt) (num_files x ord)

[Ti,g]=LM_calc_norm_str(ord,num_files,Z,Y,freq); % Συνάρτηση Levenberg-Marquardt με seeding αυξανομένης της συχνότητας που επιλύει το eigenvector switchover - Υπολογίζει πίνακα ιδιοδιανυσμάτων Ti (current) (num_files x ord^2) και ιδιοτιμές g (sqrt) (num_files x ord)

