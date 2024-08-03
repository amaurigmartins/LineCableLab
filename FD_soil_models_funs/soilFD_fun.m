function [soilFD] = soilFD_fun(soil,FD_flag,f_total)

soilFD = soil;

siz=length(f_total);

er_HF_LS=5;
er_HF_M=8;
er_HF_dry=3.5;
e0=8.854187817e-12;  % Farads/meters


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

soilFD.erg_total=erg_total;
soilFD.sigma_g_total = sigma_g_total;
soilFD.e_g_total = e_g_total;
soilFD.m_g = m_g;
end

