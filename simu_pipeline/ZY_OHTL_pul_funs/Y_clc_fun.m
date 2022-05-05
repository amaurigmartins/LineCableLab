function [Ytot_Imag,Ytot_Pet,sigma_g_total,e_g_total] = Y_clc_fun(f_total,ord,ZYprnt,FD_flag,siz,soil,h,d,Geom,ZYsave,jobid)
% Variables
e0=8.854187817e-12;  % Farads/meters
omega_total=2*pi*f_total;


%% Conductor data
% Determine the outer radius of the cable
cab_ex=zeros(ord,1);
for k=1:1:ord
    if isnan(Geom(k,8))
        cab_ex(k,1)=Geom(k,5);
    else
        cab_ex(k,1)=Geom(k,8);
    end   
end

%% Earth data
% Earth Electric Parameters
erg_total=soil.erg;         %relative permittivity of earth
m_g=soil.m_g;     %permeability of earth
sigma_g_total=soil.sigma_g; %conductivity of earth
er_HF_LS=5;
er_HF_M=8;
er_HF_dry=3.5;

rho_LF=1/sigma_g_total;

% FD soil parameters
if FD_flag == 0
    sigma_g_total=sigma_g_total.*ones(siz,1);
    erg_total=erg_total.*ones(siz,1);
elseif FD_flag == 1
    [sigma_g_total,rho_g,erg_total] = soilFD_LongSmith_fun(rho_LF,er_HF_LS,f_total);
elseif FD_flag == 2
    [sigma_g_total,rho_g,erg_total] = soilFD_Portela_fun(rho_LF,f_total); % Portela
elseif FD_flag == 3
    [sigma_g_total,rho_g,erg_total] = soilFD_AlVis_fun(rho_LF,f_total); % Alipio & Visacro
elseif FD_flag == 4
    [sigma_g_total,rho_g,erg_total] = soilFD_DM_fun(rho_LF,er_HF_dry,f_total); % Datsios & Mikropoulos
elseif FD_flag == 5    
    [sigma_g_total,rho_g,erg_total] = soilFD_Scott_fun(rho_LF,f_total); % Scott
elseif FD_flag == 6
    [sigma_g_total,rho_g,erg_total] = soilFD_M_fun(rho_LF,er_HF_M,f_total); % Messier
elseif FD_flag == 7    
    [sigma_g_total,rho_g,erg_total] = soilFD_VisPor_fun(rho_LF,f_total); % Visacro & Portela
elseif FD_flag == 8        
    [sigma_g_total,rho_g,erg_total] = soilFD_VisAl_fun(rho_LF,f_total); % Visacro & Alipio
elseif FD_flag == 9 
    [sigma_g_total,rho_g,erg_total] = soilFD_Cigre_fun(rho_LF,f_total); % Cigre 
end
e_g_total=e0.*erg_total;

%% Calculations
% Dynamic allocation
Ytot_Imag=zeros(ord,ord,siz);
Ytot_Pet=zeros(ord,ord,siz);
%% Shunt admittance
for k=1:siz
    omega=omega_total(k);
    f=f_total(k);
    erg=erg_total(k);
    sigma_g=sigma_g_total(k);
    e_g=e_g_total(k);
    n1_tetragwno=1;
%%%% Pettersson's Method %%%%

  %% Potential coefficients (self terms)

    %Potential Coeficients of Self Admittance - Perfect Ground
     Ps_pet_perf=P_pet_slf_perf(h,cab_ex,ord);
    %Potential Coeficients of Self Admittance - Imperfect Ground
     Ps_pet_imperf=P_pet_slf_imperf(h,e_g,m_g,sigma_g,omega,ord);
    % Total self - imperfect
    Ps_pet=Ps_pet_perf+Ps_pet_imperf;
    % Total self - perfect
    Ps_imag=Ps_pet_perf;
    
  %% Potential coefficients (mutual terms)
    %Potential Coeficients of Mutal Admittance - Perfect Ground
     Pm_pet_perf=P_pet_mut_perf(h,d,ord);    
    %Potential Coeficients of Mutal Admittance - Imperfect Ground
     Pm_pet_imperf=P_pet_mut_imperf(h,d,e_g,m_g,sigma_g,omega,ord);   
    % Total Mutual - imperfect
    Pm_pet=Pm_pet_perf+Pm_pet_imperf;
    % Total Mutual - imperfect
    Pm_imag=Pm_pet_perf;
 
  %%Potential Coeficients of Conductor insulation
    [Pin_mat]=Pins_mat_fun(ord,Geom);
   
    
  %% Total Pottential Coefficient Matrix
    L_Q_imag_mat=Ps_imag+Pm_imag+Pin_mat;
    L_Q_mat=Ps_pet+Pm_pet+Pin_mat;
    % Total Admittance Matrices
    Ytot_Imag(:,:,k)=1i.*omega.*e0.*2.*pi.*n1_tetragwno.*inv(L_Q_imag_mat);
    Ytot_Pet(:,:,k)=1i.*omega.*e0.*2.*pi.*n1_tetragwno.*inv(L_Q_mat);
    
end

if (ZYprnt)
    plotY_fun_ct(f_total,ord,Ytot_Imag,Ytot_Pet,jobid)
end

if (ZYsave)
fname = [jobid '_Y_param.mat'];
save(fname,'Ytot_Pet', 'Ytot_Imag', 'f');
end