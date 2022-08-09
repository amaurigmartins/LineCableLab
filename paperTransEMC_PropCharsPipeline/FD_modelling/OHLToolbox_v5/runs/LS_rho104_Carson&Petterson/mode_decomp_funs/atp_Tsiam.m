function [Ti_final,g]=atp_Tsiam(Z_tot,Y_tot,num,freq,num_files)

g=zeros(num_files,num); % (num_files x ord)
Ti_final=zeros(num_files,num^2); % (num_files x ord^2)
Zphase=zeros(num,num); % (ord x ord)
Yphase=zeros(num,num); % (ord x ord)

for a=1:1:num_files
    
    f=freq(a);
    
    for o=1:num
        Zphase(o,:)=Z_tot(a,(o-1)*num+1:o*num);
        Yphase(o,:)=Y_tot(a,(o-1)*num+1:o*num);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Modal's Algorithm Starts Here %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZginY=Zphase*Yphase;
[eig_vec,eig_val]=eig(ZginY);

for k=1:num
    eig_grammh(k)=eig_val(k,k);
end

gamma=sqrt(eig_grammh);
phase_vel=(2*pi*f)./imag(gamma);
atten=real(gamma);
ak_db=8.686.*atten;
Tv_1=eig_vec;
Y_mode_preliminary=transpose(Tv_1)*Yphase*Tv_1;
gwnia=-(angle(Y_mode_preliminary)-pi/2);
for k=1:num
    for j=1:num
        gwnia(j,k)=gwnia(k,k);
    end
end
metro=ones(num);
eak=metro.*cos(gwnia)+i.*metro.*sin(gwnia);
gwnia_dia_2=gwnia./2.0;
eak_dia_2=metro.*cos(gwnia_dia_2)+i.*metro.*sin(gwnia_dia_2);
Tv_2=eak_dia_2.*Tv_1;
Y_mode_k=transpose(Tv_2)*Yphase*Tv_2;
Tv_3=real(Tv_2);
phliko1=max(Tv_3);
phliko2=min(Tv_3);
for k=1:num
    if abs(phliko1(k))<abs(phliko2(k))
        phliko(k)=phliko2(k);
    else
        phliko(k)=phliko1(k);
    end
end
kanon=ones(num);
for k=1:num
    for j=1:num
        kanon(j,k)=phliko(k);
    end
end
Tv_4=Tv_3./kanon;
Ti_appr=inv(transpose(Tv_4));
Z_mode=inv(Tv_4)*Zphase*Ti_appr; % Possible Variable that I need!!!!
Y_mode=transpose(Tv_4)*Yphase*Tv_4; % Possible Variable that I need!!!!
for k=1:num
    R_pch(k)=real(Z_mode(k,k));
    G_pch(k)=real(Y_mode(k,k));
end
for k=1:num
    ta1(k)=Z_mode(k,k);
end
for k=1:num
    ta2(k)=Y_mode(k,k);
end

% Make Matrices for PCH or LIB Format

Zs_pch=sqrt(imag(ta1)./imag(ta2));
V_pch=phase_vel;
Tv_test=Tv_2./kanon;
Ti_pch=inv(transpose(Tv_test)); % Possible Variable that I need!!!!
Ti_pchr=real(Ti_pch);
Ti_pchi=imag(Ti_pch);
for k=1:num
    for j=1:num
        Ti(2*k-1,j)=Ti_pchr(k,j);
        Ti(2*k,j)=Ti_pchi(k,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Modal's Algorithm Ends Here %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Temporary Variables
R=R_pch;
G=G_pch;
Z=Zs_pch;
V=V_pch;
T=Ti;

Z_mode=diag(Z_mode);
Y_mode=diag(Y_mode);

real_Zmode=real(Z_mode);
imag_Zmode=imag(Z_mode);

real_Ymode=real(Y_mode);
imag_Ymode=imag(Y_mode);
Ti_real=real(Ti_pch);
Ti_imag=imag(Ti_pch);

Zchar_modes=sqrt(Z_mode./Y_mode);
gamma_modes=sqrt(Z_mode.*Y_mode);

% Characteristic Impedance Matrix in phase variables, in ohm

    
Zchar_modes_t=diag(Zchar_modes);
Zchar_phase=Tv_4*Zchar_modes_t*transpose(Tv_4);

gamma_modes=transpose(gamma_modes);

Ttemp=zeros(num,num); % (num_files x (2*ord)^2)
for k=1:num
    Ttemp(k,:)=T(2*k-1,:)+1i*T(2*k,:); % ο Ti είναι μιγαδικός!
end

g(a,:)=gamma_modes;

    for o=1:num
        Ti_final(a,(o-1)*num+1:o*num)=Ttemp(o,:);
    end
    
end
