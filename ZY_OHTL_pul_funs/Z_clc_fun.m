function [Ztot_Carson,Ztot_Noda,Ztot_Deri,Ztot_AlDe,Ztot_Sunde,Ztot_Pettersson,Ztot_Semlyen,Ztot_Wise,Ztot_Papad,Ztot_OvUnd,Ztot_Kik,Ztot_Sunde2La,Ztot_OvUndPol,Ztot_Xue,Nph] = Z_clc_fun(f_total,ord,ZYprnt,FD_flag,siz,soil,h,d,Geom,ZYsave,jobid)
%% Variables
e0=8.854187817e-12;  % Farads/meters
% m0=4*pi*1e-7;        % Henry's/meters

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
rad_in=Geom(:,4); % internal radius of conductor
rad_ex=Geom(:,5); % external radius of conductor
sigma_w=1./Geom(:,6); %conductivity of conductor
mrw=Geom(:,7); % relative permeability of conductor
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
Ztot_Carson=zeros(Nph,Nph,siz);
Ztot_Noda=zeros(Nph,Nph,siz);
Ztot_Deri=zeros(Nph,Nph,siz);
Ztot_Sunde=zeros(Nph,Nph,siz);
Ztot_Pettersson=zeros(Nph,Nph,siz);
Ztot_Wise=zeros(Nph,Nph,siz);
Ztot_Semlyen=zeros(Nph,Nph,siz);
Ztot_AlDe=zeros(Nph,Nph,siz);
Ztot_Papad=zeros(Nph,Nph,siz);
Ztot_OvUnd=zeros(Nph,Nph,siz);
Ztot_OvUndPol=zeros(Nph,Nph,siz);
Ztot_Kik=zeros(Nph,Nph,siz);
Ztot_Sunde2La=zeros(Nph,Nph,siz);
Ztot_Xue=zeros(Nph,Nph,siz);

Zin_test=zeros(ord,ord,siz);
Zg_test=zeros(ord,ord,siz);
Zg_Carson_mat=zeros(ord,ord,siz);
Zs2la=zeros(ord,ord);
Zm2la=zeros(ord,ord);

%% Series impedances
for k=1:siz
    omega=omega_total(k);
    f=f_total(k);
%     erg=erg_total(k);
    sigma_g=sigma_g_total(k);
    mu_g=m_g(end);
    e_g=e_g_total(k);
    %freq=f(k);

    Zin=Z_skin_mat_fun_ct(ord,rad_ex,rad_in,sigma_w,mrw,omega,Geom);
    Zin_test(:,:,k)=Zin;



    if all(h>0) % all conductors are aboveground

        %%%% Carson's Method %%%%
        % Formulas includes the effect of the perfect and the imperfect earth
        % Self Impedances
        Zs_carson = Z_carson_slf(h,cab_ex,sigma_g,omega,0,ord);
        % Mutual Impedances
        Zm_carson=Z_carson_mut(h,d,omega,sigma_g,0,ord);

        %%%% Noda's Method %%%%
        % Formulas includes the effect of the perfect and the imperfect earth
        % Self Impedances
        Zs_noda=Z_noda_slf(h,0,sigma_g,omega,cab_ex,ord);
        % Mutual Impedances
        Zm_noda=Z_noda_mut(h,e_g,sigma_g,omega,d,ord);

        %%%% Deri's Method %%%%
        % Self Impedances
        Zs_deri=Z_der_slf(h,mu_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_deri=Z_der_mut(h,d,mu_g,sigma_g,omega,ord);

        %%%% Alvarado - Betancourt Method %%%%
        % Self Impedances
        Zs_alde=Z_alvadetan_slf(h,e_g,mu_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_alde=Z_alvadetan_mut(h,d,e_g,mu_g,sigma_g,omega,ord);

        %%%% Sunde's Method %%%%
        % Self Impedances
        Zs_sunde=Z_snd_slf(h,e_g,mu_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_sunde=Z_snd_mut(h,d,e_g,mu_g,sigma_g,omega,ord);

        %%%% Semlyen's Method %%%%
        % Self Impedances
        Zs_semlyen=Z_sln_slf(h,e_g,mu_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_semlyen=Z_sln_mut(h,d,e_g,mu_g,sigma_g,omega,ord);

        %%%% Pettersson's Method %%%%
        % Self Impedances
        Zs_pet=Z_pet_slf(h,e_g,mu_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_pet=Z_pet_mut(h,d,e_g,mu_g,sigma_g,omega,ord);

        %%%% Wise's Integral Formula %%%%
        % Self Impedances
        Zs_wise=Z_wise_slf(h,e_g,mu_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_wise=Z_wise_mut(h,d,e_g,mu_g,sigma_g,omega,ord);

        %%%% Kikuchi's Integral Formula %%%%
        global kxa;if isempty(kxa);kxa='k0';end;
        % Self Impedances
        Zs_kik=Z_kik_slf(h,cab_ex,e_g ,mu_g,sigma_g,f,ord,kxa); % self impedances of the overhead conductors
        % Mutual Impedances
        Zm_kik=Z_kik_mut(h,d,e_g ,mu_g,sigma_g,f,ord,kxa); % mutual impedances of the overhead conductors

        if num_layers==2
            % Sunde's formula for 2-layers
            sigma_g_la=sigma_g_total(1,:);
            e_g_la=e_g_total(1,:);
            global kxa;if isempty(kxa);kxa='k0';end;
            t=-soil.layer(1).t;
            Zs2la=Z_ohl_slf_2lay(h,cab_ex,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
            Zm2la=Z_ohl_mut_2lay(h,d,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        end

        %%%% Total Matrices %%%%
        Zg_Carson=Zs_carson+Zm_carson;
        Zg_Noda=Zs_noda+Zm_noda;
        Zg_Deri=Zs_deri+Zm_deri;
        Zg_AlDe=Zs_alde+Zm_alde;
        Zg_Sund=Zs_sunde+Zm_sunde;
        Zg_Pet=Zs_pet+Zm_pet;
        Zg_Wise=Zs_wise+Zm_wise;
        Zg_Sem=Zs_semlyen+Zm_semlyen;
        Zg_kik=Zs_kik+Zm_kik;
        Zg_Sunde2La=Zs2la+Zm2la;

        Zpg=Z_self_mut_pg(h,d,cab_ex,omega,ord);        % Influence of perfect earth

        Ze_pg(:,:,k)=Zpg;
        Zg(:,:,k)=Zg_Carson-Zpg;
        Zg_Carson_mat(:,:,k)=Zg_Carson;
        Z_Carson=Zin+Zg_Carson;
        Z_Noda=Zin+Zg_Noda;
        Z_Deri=Zin+Zpg+Zg_Deri;
        Z_AlDe=Zin+Zpg+Zg_AlDe;
        Z_Sunde=Zin+Zpg+Zg_Sund;
        Z_Pettersson=Zin+Zpg+Zg_Pet;
        Z_Wise=Zin+Zpg+Zg_Wise;
        Z_kik=Zin+Zg_kik;
        Z_Sunde2La=Zin+Zg_Sunde2La;
        Z_Semlyen=Zin+Zpg+Zg_Sem;

        Zg_test(:,:,k)=Zg_Sunde2La;

        % Bundle reduction matrices
        Ztot_Carson(:,:,k) = bundleReduction(ph_order,Z_Carson);
        Ztot_AlDe(:,:,k) = bundleReduction(ph_order,Z_AlDe);
        Ztot_Deri(:,:,k) = bundleReduction(ph_order,Z_Deri);
        Ztot_Noda(:,:,k) = bundleReduction(ph_order,Z_Noda);
        Ztot_Pettersson(:,:,k) = bundleReduction(ph_order,Z_Pettersson);
        Ztot_Semlyen(:,:,k) = bundleReduction(ph_order,Z_Semlyen);
        Ztot_Sunde(:,:,k) = bundleReduction(ph_order,Z_Sunde);
        Ztot_Wise(:,:,k) = bundleReduction(ph_order,Z_Wise);
        Ztot_Kik(:,:,k) = bundleReduction(ph_order,Z_kik);
        if num_layers==2;Ztot_Sunde2La(:,:,k) = bundleReduction(ph_order,Z_Sunde2La);end

    elseif all(h<0) % all conductors are underground
        global kxe;if isempty(kxe);kxe='k1';end;
               
        Zs_papad=Z_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe);
        Zm_papad=Z_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe);
        Ztot_Papad(:,:,k)=Zs_papad+Zm_papad;

        Zs_xue=Z_xue_slf(h,cab_ex,e_g,sigma_g,f,ord);
        Zm_xue=Z_xue_mut(h,d,e_g,sigma_g,f,ord);
        Ztot_Xue(:,:,k)=Zs_xue+Zm_xue;

        Ztot_Papad(:,:,k)=Zin_test(:,:,k)+Ztot_Papad(:,:,k);
        Ztot_Xue(:,:,k)=Zin_test(:,:,k)+Ztot_Xue(:,:,k);

        Ztot_Papad(:,:,k) = bundleReduction(ph_order,Ztot_Papad(:,:,k));
        Ztot_Xue(:,:,k) = bundleReduction(ph_order,Ztot_Xue(:,:,k));
    
    else % a mixed overhead-underground setup. Need to choose carefully the self impedance formula
        
        global kxe;if isempty(kxe);kxe=0;end;
        Zs_papad=Z_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self impedances of the underground conductors
        Zm_papad=Z_papad_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxe); % mutual impedances of the underground conductors
        global kxa;if isempty(kxa);kxa=0;end;
        Zs_kik=Z_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self impedances of the overhead conductors
        Zm_kik=Z_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual impedances of the overhead conductors
        global kxm;if isempty(kxm);kxm=0;end;
        Zm_new=Z_new_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxm); % mutual impedances in the mixed configuration
        Zm_pol=Z_pol_mut(h,d,sigma_g,f,ord); % mutual impedances in the mixed configuration Pollaczek
        Ztot_OvUnd(:,:,k)=Zs_papad+Zm_papad+Zs_kik+Zm_kik+Zm_new;
        Ztot_OvUnd(:,:,k)=Zin_test(:,:,k)+Ztot_OvUnd(:,:,k);
        Ztot_OvUnd(:,:,k) = bundleReduction(ph_order,Ztot_OvUnd(:,:,k));

        Ztot_OvUndPol(:,:,k)=Zs_papad+Zm_papad+Zs_kik+Zm_kik+Zm_pol;
        Ztot_OvUndPol(:,:,k)=Zin_test(:,:,k)+Ztot_OvUndPol(:,:,k);
        Ztot_OvUndPol(:,:,k) = bundleReduction(ph_order,Ztot_OvUndPol(:,:,k));


    end


end
%% Plot parameters
if (ZYprnt)
    plotZ_fun_ct(f_total,Nph,Ztot_Carson,Ztot_Noda,Ztot_Deri,Ztot_AlDe,Ztot_Sunde,Ztot_Pettersson,Ztot_Semlyen,Ztot_Wise,Ztot_Papad,Ztot_OvUnd,Ztot_Kik,Ztot_Sunde2La,Ztot_OvUndPol,Ztot_Xue,jobid);
end

if (ZYsave)
    fname = fullfile(pwd,'Z_pul_output.mat');
    save(fname,'Ztot_Carson', 'Ztot_Noda', 'Ztot_Deri', 'Ztot_AlDe', 'Ztot_Sunde', 'Ztot_Pettersson', 'Ztot_Semlyen', 'Ztot_Wise','Ztot_Papad','Ztot_OvUnd','Ztot_Kik','Ztot_Sunde2La','Ztot_OvUndPol','Ztot_Xue');
end