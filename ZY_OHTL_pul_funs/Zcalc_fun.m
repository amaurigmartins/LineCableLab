function [out] = Zcalc_fun(f_total,ord,Nph,soilFD,h,d,Geom,opts)

if nargin == 7; opts=struct();end

% Extract all the field values from the structure
if isstruct(opts)
    fieldValues = struct2cell(opts);
else
    fieldValues={};
end

% Ensure all field values are logicals
    if ~isempty(fieldValues) && any(~cellfun(@islogical, fieldValues))
        error('All fields in the options structure must be logical values. Come on, it''s not that hard.');
    end

% Check if all fields are false
allFalse = all(~[fieldValues{:}]);

% Ensure there is at least one task to perform
if allFalse
    warning('No formula was specified. Using default.');
    if all(h>0)
        opts.Wise=true;
    elseif all(h<0)
        opts.Papadopoulos=true;
    else
        opts.OverUnder=true;
    end
end

useFormula = @(testString) isfield(opts,testString) && opts.(testString);

% Frequency samples
omega_total=2*pi*f_total;
siz=length(f_total);

% Conductor data
ph_order = Geom(:,1);
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

% Compute the actual parameters
o=1; % counter for the number of outputs
for k=1:siz
    omega=omega_total(k);
    f=f_total(k);
    sigma_g=soilFD.sigma_g_total(k);
    m_g=soilFD.m_g(end);
    e_g=soilFD.e_g_total(k);
    
    %%%% Internal impedance of the conductor
    Zin=Z_skin_mat_fun_ct(ord,rad_ex,rad_in,sigma_w,mrw,omega,Geom);

    %%%% Influence of perfect earth
    Zpg=Z_self_mut_pg(h,d,cab_ex,omega,ord);

    %%%% Carson
    if useFormula('Carson')
        if k==1; Ztot_Carson=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Formulas includes the effect of the perfect and the imperfect earth
        % Self Impedances
        Zs_carson = Z_carson_slf(h,cab_ex,sigma_g,omega,0,ord);
        % Mutual Impedances
        Zm_carson=Z_carson_mut(h,d,omega,sigma_g,0,ord);
        % Total matrices
        Zg_Carson=Zs_carson+Zm_carson;
        Z_Carson=Zin+Zg_Carson;
        Ztot_Carson(:,:,k) = bundleReduction(ph_order,Z_Carson);
        % Store outputs
        if k==siz
            out(o).VarName='Ztot_Carson';
            out(o).Label='Carson';
            out(o).Values=Ztot_Carson;
            o=o+1;
        end
    end

    %%%% Noda
    if useFormula('Noda')
        if k==1; Ztot_Noda=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Formulas includes the effect of the perfect and the imperfect earth
        % Self Impedances
        Zs_noda=Z_noda_slf(h,0,sigma_g,omega,cab_ex,ord);
        % Mutual Impedances
        Zm_noda=Z_noda_mut(h,e_g,sigma_g,omega,d,ord);
        % Total matrices
        Zg_Noda=Zs_noda+Zm_noda;
        Z_Noda=Zin+Zg_Noda;
        Ztot_Noda(:,:,k) = bundleReduction(ph_order,Z_Noda);
        % Store outputs
        if k==siz
            out(o).VarName='Ztot_Noda';
            out(o).Label='Noda';
            out(o).Values=Ztot_Noda;
            o=o+1;
        end
    end
    
    %%%% Deri
    if useFormula('Deri')
        if k==1; Ztot_Deri=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self Impedances
        Zs_deri=Z_der_slf(h,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_deri=Z_der_mut(h,d,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_Deri=Zpg+Zs_deri+Zm_deri;
        Z_Deri=Zin+Zg_Deri;
        Ztot_Deri(:,:,k) = bundleReduction(ph_order,Z_Deri);
        if k==siz
            out(o).VarName='Ztot_Deri';
            out(o).Label='Deri';
            out(o).Values=Ztot_Deri;
            o=o+1;
        end
    end

    %%%% Alvarado - Betancourt
    if useFormula('Alvarado')
        if k==1; Ztot_AlBe=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self Impedances
        Zs_alde=Z_alvadetan_slf(h,e_g,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_alde=Z_alvadetan_mut(h,d,e_g,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_AlBe=Zpg+Zs_alde+Zm_alde;
        Z_AlBe=Zin+Zg_AlBe;
        Ztot_AlBe(:,:,k) = bundleReduction(ph_order,Z_AlBe);
        if k==siz
            out(o).VarName='Ztot_AlBe';
            out(o).Label='Alvarado-Betancourt';
            out(o).Values=Ztot_AlBe;
            o=o+1;
        end
    end
    
    %%%% Sunde (uniform soil)
    if useFormula('Sunde')
        if k==1; Ztot_Sunde=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self Impedances
        Zs_sunde=Z_snd_slf(h,e_g,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_sunde=Z_snd_mut(h,d,e_g,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_Sund=Zs_sunde+Zm_sunde+Zpg;
        Z_Sunde=Zin+Zg_Sund;
        Ztot_Sunde(:,:,k) = bundleReduction(ph_order,Z_Sunde);
        if k==siz
            out(o).VarName='Ztot_Sunde';
            out(o).Label='Sunde (uniform)';
            out(o).Values=Ztot_Sunde;
            o=o+1;
        end
    end
    
    %%%% Semlyen
    if useFormula('Semlyen')
        if k==1; Ztot_Semlyen=zeros(Nph,Nph,siz); end % Prelocate matrix
         % Self Impedances
        Zs_semlyen=Z_sln_slf(h,e_g,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_semlyen=Z_sln_mut(h,d,e_g,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_Sem=Zs_semlyen+Zm_semlyen+Zpg;
        Z_Semlyen=Zin+Zg_Sem;
        Ztot_Semlyen(:,:,k) = bundleReduction(ph_order,Z_Semlyen);
        if k==siz
            out(o).VarName='Ztot_Semlyen';
            out(o).Label='Semlyen';
            out(o).Values=Ztot_Semlyen;
            o=o+1;
        end
    end
    
    %%%% Pettersson
    if useFormula('Pettersson')
        if k==1; Ztot_Pettersson=zeros(Nph,Nph,siz); end % Prelocate matrix
         % Self Impedances
        Zs_pet=Z_pet_slf(h,e_g,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_pet=Z_pet_mut(h,d,e_g,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_Pet=Zs_pet+Zm_pet+Zpg;
        Z_Pettersson=Zin+Zg_Pet;
        Ztot_Pettersson(:,:,k) = bundleReduction(ph_order,Z_Pettersson);
        if k==siz
            out(o).VarName='Ztot_Pettersson';
            out(o).Label='Pettersson';
            out(o).Values=Ztot_Pettersson;
            o=o+1;
        end
    end

    %%%% Wise
    if useFormula('Wise')
        if k==1; Ztot_Wise=zeros(Nph,Nph,siz); end % Prelocate matrix
         % Self Impedances
        Zs_wise=Z_wise_slf(h,e_g,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_wise=Z_wise_mut(h,d,e_g,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_Wise=Zs_wise+Zm_wise+Zpg;
        Z_Wise=Zin+Zg_Wise;
        Ztot_Wise(:,:,k) = bundleReduction(ph_order,Z_Wise);
        if k==siz
            out(o).VarName='Ztot_Wise';
            out(o).Label='Wise';
            out(o).Values=Ztot_Wise;
            o=o+1;
        end
    end

    %%%% Kikuchi
    if useFormula('Kikuchi')
        if k==1; Ztot_Kik=zeros(Nph,Nph,siz); end % Prelocate matrix
        global kxa;if isempty(kxa);kxa='k0';end;
         % Self Impedances
        Zs_kik=Z_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self impedances of the overhead conductors
        % Mutual Impedances
        Zm_kik=Z_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual impedances of the overhead conductors
        % Total matrices
        Zg_kik=Zs_kik+Zm_kik;
        Z_kik=Zin+Zg_kik;
        Ztot_Kik(:,:,k) = bundleReduction(ph_order,Z_kik);
        if k==siz
            out(o).VarName='Ztot_Kik';
            out(o).Label='Kikuchi';
            out(o).Values=Ztot_Kik;
            o=o+1;
        end
    end

    %%%% Sunde (2-layered soil)
    if useFormula('Sunde2Layers')
        if k==1; Ztot_Sunde2La=zeros(Nph,Nph,siz); end % Prelocate matrix
        sigma_g_la=soilFD.sigma_g_total(1,:);
        e_g_la=soilFD.e_g_total(1,:);
        global kxa;if isempty(kxa);kxa='k0';end;
        t=-soilFD.layer(1).t;
        Zs2la=Z_ohl_slf_2lay(h,cab_ex,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        Zm2la=Z_ohl_mut_2lay(h,d,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        % Total matrices
        Zg_Sunde2La=Zs2la+Zm2la;
        Z_Sunde2La=Zin+Zg_Sunde2La;
        Ztot_Sunde2La(:,:,k) = bundleReduction(ph_order,Z_Sunde2La);
        if k==siz
            out(o).VarName='Ztot_Sunde2La';
            out(o).Label='Sunde (2-layers)';
            out(o).Values=Ztot_Sunde2La;
            o=o+1;
        end
    end

    %%%% Papadopoulos (underground)
    if useFormula('Papadopoulos')
        if k==1; Ztot_Papad=zeros(Nph,Nph,siz); end % Prelocate matrix
        global kxe;if isempty(kxe);kxe='k1';end;
        Zs_papad=Z_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe);
        Zm_papad=Z_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe);
        % Total matrices
        Zg_Papad=Zs_papad+Zm_papad;
        Z_Papad=Zin+Zg_Papad;
        Ztot_Papad(:,:,k) = bundleReduction(ph_order,Z_Papad);
        if k==siz
            out(o).VarName='Ztot_Papad';
            out(o).Label='Papadopoulos (underground)';
            out(o).Values=Ztot_Papad;
            o=o+1;
        end
    end
    
    %%%% Pollaczek (underground)
    if useFormula('Pollaczek')
        if k==1; Ztot_Pol=zeros(Nph,Nph,siz); end % Prelocate matrix
        global kxe;if isempty(kxe);kxe=0;end;
        Zs_pol=Z_papad_slf(h,cab_ex,0*e_g ,m_g,sigma_g,f,ord,kxe); % Pollackzek's result is attained by setting e_g=0
        Zm_pol=Z_papad_mut(h,d,0*e_g,m_g,sigma_g,f,ord,kxe); % Pollackzek's result is attained by setting e_g=0
        % Total matrices
        Zg_pol=Zs_pol+Zm_pol;
        Z_Pol=Zin+Zg_pol;
        Ztot_Pol(:,:,k) = bundleReduction(ph_order,Z_Pol);
        if k==siz
            out(o).VarName='Ztot_Pol';
            out(o).Label='Pollaczek (underground)';
            out(o).Values=Ztot_Papad;
            o=o+1;
        end
    end
    
    %%%% Xue (underground)
    if useFormula('Xue')
        if k==1; Ztot_Xue=zeros(Nph,Nph,siz); end % Prelocate matrix
        Zs_xue=Z_xue_slf(h,cab_ex,e_g,sigma_g,f,ord);
        Zm_xue=Z_xue_mut(h,d,e_g,sigma_g,f,ord);
        % Total matrices
        Zg_xue=Zs_xue+Zm_xue;
        Z_Xue=Zin+Zg_xue;
        Ztot_Xue(:,:,k) = bundleReduction(ph_order,Z_Xue);
        if k==siz
            out(o).VarName='Ztot_Xue';
            out(o).Label='Xue (underground)';
            out(o).Values=Ztot_Xue;
            o=o+1;
        end
    end
    
    %%%% Martins-Papadopoulos-Chrysochos (overhead-underground)
    if useFormula('OverUnder')
        if k==1; Ztot_OverUnder=zeros(Nph,Nph,siz); end % Prelocate matrix
        global kxe;if isempty(kxe);kxe=0;end;
        Zs_papad=Z_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self impedances of the underground conductors
        Zm_papad=Z_papad_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxe); % mutual impedances of the underground conductors
        global kxa;if isempty(kxa);kxa=0;end;
        Zs_kik=Z_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self impedances of the overhead conductors
        Zm_kik=Z_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual impedances of the overhead conductors
        global kxm;if isempty(kxm);kxm=0;end;
        Zm_new=Z_new_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxm); % mutual impedances in the mixed configuration
        % Total matrices
        Zg_OverUnder=Zs_papad+Zm_papad+Zs_kik+Zm_kik+Zm_new;
        Z_OverUnder=Zin+Zg_OverUnder;
        Ztot_OverUnder(:,:,k) = bundleReduction(ph_order,Z_OverUnder);
        if k==siz
            out(o).VarName='Ztot_OverUnder';
            out(o).Label='Martins-Papadopoulos-Chrysochos (overhead-underground)';
            out(o).Values=Ztot_OverUnder;
            o=o+1;
        end
    end
    
    %%%% Carson-Pollaczek (overhead-underground)
    if useFormula('CarsonPol')
        if k==1; Ztot_CarsonPol=zeros(Nph,Nph,siz); end % Prelocate matrix
        Zm_pol_ovund=Z_pol_mut(h,d,sigma_g,f,ord); % mutual impedances in the mixed configuration Pollaczek
        Zs_pol_und=Z_papad_slf(h,cab_ex,0*e_g ,m_g,sigma_g,f,ord,kxe); % self impedances of the underground conductors - Pollackzek's result is attained by setting e_g=0
        Zm_pol_und=Z_papad_mut(h,d,e_g ,0*m_g,sigma_g,f,ord,kxe); % mutual impedances of the underground conductors
        Zs_pol_ov=Z_kik_slf(h,cab_ex,0*e_g ,m_g,sigma_g,f,ord,kxa); % self impedances of the overhead conductors - Carson's result is attained by setting e_g=0
        Zm_pol_ov=Z_kik_mut(h,d,0*e_g ,m_g,sigma_g,f,ord,kxa); % mutual impedances of the overhead conductors
        % Total matrices
        Zg_CarsonPol=Zs_pol_und+Zm_pol_und+Zs_pol_ov+Zm_pol_ov+Zm_pol_ovund;
        Z_CarsonPol=Zin+Zg_CarsonPol;
        Ztot_CarsonPol(:,:,k) = bundleReduction(ph_order,Z_CarsonPol);
        if k==siz
            out(o).VarName='Ztot_CarsonPol';
            out(o).Label='Carson-Pollaczek (overhead-underground)';
            out(o).Values=Ztot_CarsonPol;
            o=o+1;
        end
    end
end
% end of main loop

end

