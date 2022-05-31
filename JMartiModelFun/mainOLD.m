close all
clear
clc

addpath('vfit3')

% Load dataset
fname='carson_cp_ohlt_full';
load([fname '.mat'])

% General parameters
Ns = length(f); % Number of samples
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency
line_length = 1000;
f_Ti=200e3; % Frequency chosen for JMarti model

% Raw impedance and admittance matrices
Z=permute(eval(Zmod_src),[3 1 2]);
Y=permute(eval(Ymod_src),[3 1 2]);

% Transform matrix from OHLT
clear T_temp g_temp Ti g_dis
[Ti,g_dis]=LM_calc_norm_str(ord,freq_siz,Z,Y,f);
Ti=rot_min_imag(Ti,ord,freq_siz);
for k=1:1:freq_siz
    for o=1:ord
        T_temp(o,:,k)=Ti(k,(o-1)*ord+1:o*ord);
    end
end
T.var=permute(T_temp,[3 1 2]);
g.var=g_dis;

% original transform matrix from OHTL @ chosen frequency
clear T_temp g_temp Ti g_dis
T_temp=squeeze(interp1(f,T.var,f_Ti)); % Interpolate matrix at the chosen freq
T_temp=real(T_temp); %take only the real part because of reasons...
T_temp=repmat(T_temp,1,1,freq_siz); %repeat so we can use the same original codes from OHTL
T_temp=permute(T_temp,[3 1 2]); %make frequency the first dimension
T.lm=T_temp;
g_temp=interp1(f,g.var,f_Ti);
g_temp=squeeze(repmat(g_temp,1,1,freq_siz)); 
g.lm=g_temp.';

% modified transform matrix using regular eig() from matlab
clear T_temp g_temp Ti g_dis l_temp
YY=squeeze(interp1(f,Y,f_Ti));
ZZ=squeeze(interp1(f,Z,f_Ti)); 
[T_temp,l_temp]=eig(YY*ZZ);
T_temp=real(rot(T_temp));
T_temp=repmat(T_temp,1,1,freq_siz); %repeat so we can use the same original codes from OHTL
T_temp=permute(T_temp,[3 1 2]); %make frequency the first dimension
T.eig=T_temp;
g_temp=sqrt(diag(l_temp)).';
g_temp=squeeze(repmat(g_temp,1,1,freq_siz)); 
g.eig=g_temp.';

% yet another modified version using eig() and  eigenvalue tracing
clear T_temp g_temp Ti g_dis
[Ti,g_dis]=intercheig_QR_decomp(ord,freq_siz,Z,Y);
Ti=rot_min_imag(Ti,ord,freq_siz);
for k=1:1:freq_siz
    for o=1:ord
        T_temp(o,:,k)=Ti(k,(o-1)*ord+1:o*ord);
    end
end
T_temp=permute(T_temp,[3 1 2]);

T_temp=squeeze(interp1(f,T_temp,f_Ti)); % Interpolate matrix at the chosen freq
T_temp=real(T_temp); %take only the real part because of reasons...
T_temp=repmat(T_temp,1,1,freq_siz); %repeat so we can use the same original codes from OHTL
T_temp=permute(T_temp,[3 1 2]); %make frequency the first dimension
T.intercheig=T_temp;
g_temp=interp1(f,g_dis,f_Ti);
g_temp=squeeze(repmat(g_temp,1,1,freq_siz)); 
g.intercheig=g_temp.';

clc
display('Levenberg')
squeeze(T.lm(end,:,:))
display('eig')
squeeze(T.eig(end,:,:))
display('intercheig')
squeeze(T.intercheig(end,:,:))
display('huh?')
squeeze(real(T.var(end,:,:)))

% Recompute modal params using modified T
Tsrc=T.var;
gsrc=g.var;
[Zch_m,Ych_m,~,~]=calc_char_imped_admit(Tsrc,Z,Y,ord,freq_siz);
[Zch_mod,Ych_mod,~,~]=calc_char_imped_admit(T.var,Z,Y,ord,freq_siz);

% Find optimum delay time
[tau0, tau_mps, tau_opt]=findoptimtau(f,gsrc,line_length,ord);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%DEBUGME
figure
m=6;
semilogx(f,abs(Zch_m(:,m)));hold all; semilogx(f,abs(Zch_mod(:,m)),'o')
legend('Const T @ 200 kHz (LM)', 'Var Ti')
title(sprintf('Zc mode #%d',m))



