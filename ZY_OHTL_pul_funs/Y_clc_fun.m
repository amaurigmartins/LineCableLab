function [Ytot_Imag,Ytot_Pet,Ytot_Wise,Ytot_Papad,Ytot_OvUnd,Ytot_Kik,Ytot_Sunde2La,Ytot_Xue,sigma_g_total,erg_total,Nph] = Y_clc_fun(f_total,ord,ZYprnt,FD_flag,siz,soil,h,d,Geom,ZYsave,jobid)
% Variables
e0=8.854187817e-12;  % Farads/meters
omega_total=2*pi*f_total;


ph_order = Geom(:,1);
Nph = unique(ph_order);
Nph = size(Nph(Nph~=0),1);

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
er_HF_LS=5;
er_HF_M=8;
er_HF_dry=3.5;

num_layers=1;
if isfield(soil,'layer')
    num_layers=num_layers+length(soil.layer);
end

for l=1:num_layers
    % Earth Electric Parameters
    if l==num_layers
        epsr(:,l)=soil.erg;         %relative permittivity of earth
        m_g(:,l)=soil.m_g;     %permeability of earth
        sig(:,l)=soil.sigma_g; %conductivity of earth
    elseif l < num_layers
        epsr(:,l)=soil.layer(l).erg;         %relative permittivity of earth
        m_g(:,l)=soil.layer(l).m_g;     %permeability of earth
        sig(:,l)=soil.layer(l).sigma_g; %conductivity of earth
    end

    rho_LF(:,l)=1./sig(:,l);

    % FD soil parameters
    if FD_flag == 0
        sigma_g_total(:,l)=sig(:,l).*ones(siz,1);
        erg_total(:,l)=epsr(:,l).*ones(siz,1);
    elseif FD_flag == 1
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_LongSmith_fun(rho_LF(:,l),er_HF_LS,f_total);
    elseif FD_flag == 2
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_Portela_fun(rho_LF(:,l),f_total); % Portela
    elseif FD_flag == 3
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_AlVis_fun(rho_LF(:,l),f_total); % Alipio & Visacro
    elseif FD_flag == 4
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_DM_fun(rho_LF(:,l),er_HF_dry,f_total); % Datsios & Mikropoulos
    elseif FD_flag == 5
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_Scott_fun(rho_LF(:,l),f_total); % Scott
    elseif FD_flag == 6
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_M_fun(rho_LF(:,l),er_HF_M,f_total); % Messier
    elseif FD_flag == 7
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_VisPor_fun(rho_LF(:,l),f_total); % Visacro & Portela
    elseif FD_flag == 8
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_VisAl_fun(rho_LF(:,l),f_total); % Visacro & Alipio
    elseif FD_flag == 9
        [sigma_g_total(:,l),rho_g(:,l),erg_total(:,l)] = soilFD_Cigre_fun(rho_LF(:,l),f_total); % Cigre
    end
    e_g_total(:,l)=e0.*erg_total(:,l);

end

%% Calculations
% Dynamic allocation
Ytot_Imag=zeros(Nph,Nph,siz);
Ytot_Pet=zeros(Nph,Nph,siz);
Ytot_Wise=zeros(Nph,Nph,siz);
Ptot_Papad=zeros(Nph,Nph,siz);
Ptot_OvUnd=zeros(Nph,Nph,siz);
Ytot_Papad=zeros(Nph,Nph,siz);
Ytot_OvUnd=zeros(Nph,Nph,siz);
Ytot_Kik=zeros(Nph,Nph,siz);
Ytot_Sunde2La=zeros(Nph,Nph,siz);
Ytot_Xue=zeros(Nph,Nph,siz);
Ptot_Xue=zeros(Nph,Nph,siz);
Pin_test=zeros(ord,ord,siz);
Pg_test=zeros(ord,ord,siz);
Ps2la=zeros(ord,ord);
Pm2la=zeros(ord,ord);

%% Shunt admittance
for k=1:siz
    omega=omega_total(k);
    f=f_total(k);
%     erg=erg_total(k);
    sigma_g=sigma_g_total(k);
    mu_g=m_g(end);
    e_g=e_g_total(k);
    n1_tetragwno=1;
    
        %%Potential Coeficients of Conductor insulation
    [Pin_mat]=Pins_mat_fun(ord,Geom);
     Pin_test(:,:,k)=Pin_mat./(2*pi*e0);


    if all(h>0) % all conductors are aboveground

    % Pettersson's Method
    %*** Potential coefficients (self terms)
    %Potential Coeficients of Self Admittance - Perfect Ground
    Ps_pet_perf=P_pet_slf_perf(h,cab_ex,ord);
    %Potential Coeficients of Self Admittance - Imperfect Ground
    Ps_pet_imperf=P_pet_slf_imperf(h,e_g,mu_g,sigma_g,omega,ord);
    % Total self - imperfect
    Ps_pet=Ps_pet_perf+Ps_pet_imperf;
    % Total self - perfect
    Ps_imag=Ps_pet_perf;

    % Wise's formula
    Ps_wise_imperf=P_wise_slf_imperf(h,e_g,mu_g,sigma_g,omega,ord);
    Ps_wise=Ps_pet_perf+Ps_wise_imperf;

    

    %% Potential coefficients (mutual terms)
    %Potential Coeficients of Mutal Admittance - Perfect Ground
    Pm_pet_perf=P_pet_mut_perf(h,d,ord);
    %Potential Coeficients of Mutal Admittance - Imperfect Ground
    Pm_pet_imperf=P_pet_mut_imperf(h,d,e_g,mu_g,sigma_g,omega,ord);
    % Total Mutual - imperfect
    Pm_pet=Pm_pet_perf+Pm_pet_imperf;
    % Total Mutual - imperfect
    Pm_imag=Pm_pet_perf;
    
    % Wise's formula
    Pm_wise_imperf=P_wise_mut_imperf(h,d,e_g,mu_g,sigma_g,omega,ord);
    Pm_wise=Pm_pet_perf+Pm_wise_imperf;

    % Kikuchi's formula - Self and mutual
    global kxa;if isempty(kxa);kxa='k0';end;
    Ps_kik=P_kik_slf(h,cab_ex,e_g ,mu_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
    Pm_kik=P_kik_mut(h,d,e_g,mu_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors

    if num_layers==2
        % Sunde's formula for 2-layers
        sigma_g_la=sigma_g_total(1,:);
        e_g_la=e_g_total(1,:);
        global kxa;if isempty(kxa);kxa='k0';end;
        t=-soil.layer(1).t;
        Ps2la=P_ohl_slf_2lay(h,cab_ex,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        Pm2la=P_ohl_mut_2lay(h,d,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        Pg_test(:,:,k)=Ps2la+Pm2la;
    end
   
    
    %% Total Pottential Coefficient Matrix
    L_Q_imag_mat=Ps_imag+Pm_imag+Pin_mat;
    L_Q_mat=Ps_pet+Pm_pet+Pin_mat;
    L_Q_wise=Ps_wise+Pm_wise+Pin_mat;
    L_Q_kik=Ps_kik+Pm_kik+(Pin_mat./(2*pi*e0));
    if num_layers==2;L_Q_Sunde2La=Ps2la+Pm2la+(Pin_mat./(2*pi*e0));end
    
    % Bundle reduction
    L_Q_imag_mat = bundleReduction(ph_order,L_Q_imag_mat);
    L_Q_mat = bundleReduction(ph_order,L_Q_mat);
    L_Q_wise = bundleReduction(ph_order,L_Q_wise);
    L_Q_kik = bundleReduction(ph_order,L_Q_kik);
    if num_layers==2;L_Q_Sunde2La = bundleReduction(ph_order,L_Q_Sunde2La);end
    
    % Total Admittance Matrices
    Ytot_Imag(:,:,k)=1i.*omega.*e0.*2.*pi.*n1_tetragwno.*inv(L_Q_imag_mat);
    Ytot_Pet(:,:,k)=1i.*omega.*e0.*2.*pi.*n1_tetragwno.*inv(L_Q_mat);
    Ytot_Wise(:,:,k)=1i.*omega.*e0.*2.*pi.*n1_tetragwno.*inv(L_Q_wise);
    Ytot_Kik(:,:,k)=1i.*omega.*inv(L_Q_kik);
    if num_layers==2;Ytot_Sunde2La(:,:,k)=1i.*omega.*inv(L_Q_Sunde2La);end

    elseif all(h<0) % all conductors are underground
        global kxe;if isempty(kxe);kxe='k1';end;
        Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe);
        Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe);
        Ptot_Papad(:,:,k)=Ps_papad+Pm_papad+(Pin_mat./(2*pi*e0));
        Ptot_Papad(:,:,k) = bundleReduction(ph_order,Ptot_Papad(:,:,k));

        Ps_xue=P_xue_slf(h,cab_ex,e_g,sigma_g,f,ord);
        Pm_xue=P_xue_mut(h,d,e_g,sigma_g,f,ord);
        Ptot_Xue(:,:,k)=Ps_xue+Pm_xue+(Pin_mat./(2*pi*e0));
        Ptot_Xue(:,:,k) = bundleReduction(ph_order,Ptot_Xue(:,:,k));

        Ytot_Papad(:,:,k)=1i.*omega.*inv(Ptot_Papad(:,:,k));
        Ytot_Xue(:,:,k)=1i.*omega.*inv(Ptot_Xue(:,:,k));
    
    else %over-under
        global kxe;if isempty(kxe);kxe=0;end;
        %set 'k0' (string) to air layer, 'k1' for earth, 0 to neglect
        Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self coefficients of the underground conductors
        Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe); % mutual coefficients of the underground conductors
        global kxa;if isempty(kxa);kxa=0;end;
        Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
        Pm_kik=P_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        global kxm;if isempty(kxm);kxm=0;end;
        global Pm_flag;if isempty(Pm_flag);Pm_flag=true;end;
        if Pm_flag
            Pm_new=P_new_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxm); % mutual coefficients in the mixed configuration
        else
            Pm_new=zeros(ord,ord);
        end
        Ptot_OvUnd(:,:,k)=Ps_papad+Pm_papad+Ps_kik+Pm_kik+Pm_new+(Pin_mat./(2*pi*e0));
        Ptot_OvUnd(:,:,k) = bundleReduction(ph_order,Ptot_OvUnd(:,:,k));
        Ytot_OvUnd(:,:,k)=1i.*omega.*inv(Ptot_OvUnd(:,:,k));
        
        Ptot_OvUndClassicalTL(:,:,k)=Ps_papad+Pm_papad+Ps_kik+Pm_kik+(Pin_mat./(2*pi*e0));
        Ptot_OvUndClassicalTL(:,:,k) = bundleReduction(ph_order,Ptot_OvUndClassicalTL(:,:,k));
        Ytot_OvUndClassicalTL(:,:,k)=1i.*omega.*inv(Ptot_OvUndClassicalTL(:,:,k));

    end
end

% Plot stuff
Y_pul.Ytot_Imag=Ytot_Imag;
Y_pul.Ytot_Pet=Ytot_Pet;
Y_pul.Ytot_Wise=Ytot_Wise;
Y_pul.Ytot_Papad=Ytot_Papad;
Y_pul.Ytot_OvUnd=Ytot_OvUnd;
Y_pul.Ytot_Kik=Ytot_Kik;
Y_pul.Ytot_Sunde2La=Ytot_Sunde2La;
Y_pul.Ytot_Xue=Ytot_Xue;
Y_pul.Ytot_OvUndClassicalTL=Ytot_OvUndClassicalTL;


if (ZYprnt)
    plotY_fun_ct(f_total,Nph,Y_pul,jobid)
end

if (ZYsave)
    fname = fullfile(pwd, 'Y_pul_output.mat');
    save(fname,'Ytot_Pet', 'Ytot_Imag', 'Ytot_Wise', 'Ytot_Papad', 'Ytot_OvUnd', 'Ytot_Kik','Ytot_Sunde2La','Ytot_Xue');
end