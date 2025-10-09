function [out] = Zcalc_fun(f_total,ord,Nph,soilFD,h,d,Geom,currPath,jobid,opts)

if nargin == 9; opts=struct();end

% Extract all the field values from the structure
fieldNames=fieldnames(opts);
idx=find(~(strcmp(fieldNames, 'CYZfile') | strcmp(fieldNames, 'MATfile') | strcmp(fieldNames, 'EHEMflag')));

islayered_earth = soilFD.num_layers > 1;
is_overhead = all(h>0);
is_underground = all(h<0);

if isstruct(opts)
    fieldValues = struct2cell(opts);
    fieldValues = fieldValues(idx);
else
    fieldValues={};
end

% Ensure all field values are logicals
if ~isempty(fieldValues) && any(~cellfun(@islogical, fieldValues))
    error('All fields in the formula options must be logical values. Come on, it''s not that hard.');
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
% erw=Geom(:,10); % relative permittivity of insulation
% dc=cab_ex-rad_ex; % coating thickness

% Compute the actual parameters
o=1; % counter for the number of outputs
m_g=soilFD.m_g;
for k=1:siz
    omega=omega_total(k);
    f=f_total(k);
    % uniform soil formulations use sigma_g and e_g of the bottom layer by default
    sigma_g=soilFD.sigma_g_total(k,end);
    e_g=soilFD.e_g_total(k,end);
    append_tag = 'bottom layer';
    % uses EHEM approach for layered soils
    if soilFD.num_layers > 1
        sigma_g_la=soilFD.sigma_g_total(k,:);
        e_g_la=soilFD.e_g_total(k,:);
        m_g_la=[soilFD.layer(:).m_g soilFD.m_g];
        if isfield(opts,'EHEMflag')
            if opts.EHEMflag == 1 && isfield(soilFD,'sigma_eff')
                sigma_g=soilFD.sigma_eff(k);
                append_tag = 'EHEM-sigma';
            elseif opts.EHEMflag == 2 && isfield(soilFD,'sigma_eff') && isfield(soilFD,'eps_eff')
                sigma_g=soilFD.sigma_eff(k);
                e_g=soilFD.eps_eff(k);
                append_tag = 'EHEM-gamma';
            elseif opts.EHEMflag == 11 && isfield(soilFD,'sigma_eff') % artificial homogeneous earth implemented via identical layers
                sigma_g=soilFD.sigma_eff(k);
                sigma_g_la(1,1:soilFD.num_layers)=repmat(soilFD.sigma_eff(k),1,soilFD.num_layers);
                e_g_la(1,1:soilFD.num_layers)=repmat(e_g,1,soilFD.num_layers);
                append_tag = 'EHEM-sigma';
            elseif opts.EHEMflag == 21 && isfield(soilFD,'sigma_eff') && isfield(soilFD,'eps_eff') % artificial homogeneous earth implemented via identical layers
                sigma_g=soilFD.sigma_eff(k);
                e_g=soilFD.eps_eff(k);
                sigma_g_la(1,1:soilFD.num_layers)=repmat(soilFD.sigma_eff(k),1,soilFD.num_layers);
                e_g_la(1,1:soilFD.num_layers)=repmat(soilFD.eps_eff(k),1,soilFD.num_layers);
                append_tag = 'EHEM-gamma';
            end
        end
    end

    %%%% Internal impedance of the conductor
    Zin=Z_skin_mat_fun_ct(ord,rad_ex,rad_in,sigma_w,mrw,omega,Geom);

    %%%% Influence of perfect earth
    Zpg=Z_self_mut_pg(h,d,cab_ex,omega,ord);

    %%%% Carson
    if useFormula('Carson') && is_overhead
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
    if useFormula('Noda') && is_overhead
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
    if useFormula('Deri') && is_overhead
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
    if useFormula('Alvarado') && is_overhead
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
    if useFormula('Sunde') && is_overhead
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
            out(o).Label='Sunde';
            out(o).Values=Ztot_Sunde;
            o=o+1;
        end
    end

    %%%% Semlyen
    if useFormula('Semlyen') && is_overhead
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
    if useFormula('Pettersson') && is_overhead
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
    if useFormula('Wise') && is_overhead
        if k==1; Ztot_Wise=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self Impedances
        Zs_wise=Z_wise_slf(h,cab_ex,e_g,m_g,sigma_g,omega,ord);
        % Mutual Impedances
        Zm_wise=Z_wise_mut(h,d,e_g,m_g,sigma_g,omega,ord);
        % Total matrices
        Zg_Wise=Zs_wise+Zm_wise; %+Zpg;
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
    % if useFormula('Kikuchi') && is_overhead
    %     if k==1; Ztot_Kik=zeros(Nph,Nph,siz); end % Prelocate matrix
    %     kxa='k0';
    %     % Self Impedances
    %     Zs_kik=Z_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self impedances of the overhead conductors
    %     % Mutual Impedances
    %     Zm_kik=Z_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual impedances of the overhead conductors
    %     % Total matrices
    %     Zg_kik=Zs_kik+Zm_kik;
    %     Z_kik=Zin+Zg_kik;
    %     Ztot_Kik(:,:,k) = bundleReduction(ph_order,Z_kik);
    %     if k==siz
    %         out(o).VarName='Ztot_Kik';
    %         out(o).Label='Kikuchi';
    %         out(o).Values=Ztot_Kik;
    %         o=o+1;
    %     end
    % end

    %%%% Nakagawa (2-layered soil)
    if useFormula('Naka2Layers') && islayered_earth && is_overhead
        if k==1; Ztot_Naka2La=zeros(Nph,Nph,siz); end % Prelocate matrix
        % sigma_g_la=soilFD.sigma_g_total(1,:);
        % e_g_la=soilFD.e_g_total(1,:);
        % m_g_la=[soilFD.layer(:).m_g soilFD.m_g];
        %0 for sunde, k0 for nakagawa %%%%%%%%%%%%%%%%%%% CHECK ME
        kxa='k0'; 
        t=-soilFD.layer(1).t;
        Zs2la=Z_ohl_slf_2lay(h,cab_ex,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxa);
        Zm2la=Z_ohl_mut_2lay(h,d,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxa);
        % Total matrices
        Zg_Naka2La=Zs2la+Zm2la;
        Z_Naka2La=Zin+Zg_Naka2La;
        Ztot_Naka2La(:,:,k) = bundleReduction(ph_order,Z_Naka2La);
        if k==siz
            out(o).VarName='Ztot_Naka2La';
            out(o).Label='Nakagawa (overhead, 2-layers)';
            out(o).Values=Ztot_Naka2La;
            o=o+1;
        end
    end

    %%%% Papadopoulos (underground, 2-layered soil)
    if useFormula('Papad2LayersUnder') && islayered_earth && is_underground
        if k==1; Ztot_Papad2LaUnder=zeros(Nph,Nph,siz); end % Prelocate matrix
        % sigma_g_la=soilFD.sigma_g_total(1,:);
        % e_g_la=soilFD.e_g_total(1,:);
        % m_g_la=[soilFD.layer(:).m_g soilFD.m_g];
        kxe=0;
        t=-soilFD.layer(1).t;
        Zs2la_under=Z_papad_slf_2lay_under(h,cab_ex,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxe);
        Zm2la_under=Z_papad_mut_2lay_under(h,d,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxe);
        % Total matrices
        Zg_Papad2LaUnder=Zs2la_under+Zm2la_under;
        Z_Papad2LaUnder=Zin+Zg_Papad2LaUnder;
        Ztot_Papad2LaUnder(:,:,k) = bundleReduction(ph_order,Z_Papad2LaUnder);
        if k==siz
            out(o).VarName='Ztot_Papad2LaUnder';
            out(o).Label='Papadopoulos (underground, 2-layers)';
            out(o).Values=Ztot_Papad2LaUnder;
            o=o+1;
        end
    end
   
    %%%% Papadopoulos (underground)
    if useFormula('Papadopoulos') && is_underground
        if k==1; Ztot_Papad=zeros(Nph,Nph,siz); end % Prelocate matrix
        kxe='k1';
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
    if useFormula('Pollaczek') && is_underground
        if k==1; Ztot_Pol=zeros(Nph,Nph,siz); end % Prelocate matrix
        kxe=0;
        Zs_pol=Z_papad_slf(h,cab_ex,0*e_g,m_g,sigma_g,f,ord,kxe); % Pollackzek's result is attained by setting e_g=0
        Zm_pol=Z_papad_mut(h,d,0*e_g,m_g,sigma_g,f,ord,kxe); % Pollackzek's result is attained by setting e_g=0
        % Total matrices
        Zg_pol=Zs_pol+Zm_pol;
        Z_Pol=Zin+Zg_pol;
        Ztot_Pol(:,:,k) = bundleReduction(ph_order,Z_Pol);
        if k==siz
            out(o).VarName='Ztot_Pol';
            out(o).Label='Pollaczek (underground)';
            out(o).Values=Ztot_Pol;
            o=o+1;
        end
    end

    %%%% Xue (underground)
    if useFormula('Xue') && is_underground
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
        kxe=0;
        Zs_papad=Z_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self impedances of the underground conductors
        Zm_papad=Z_papad_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxe); % mutual impedances of the underground conductors
        Zs_over=Z_wise_slf(h,cab_ex,e_g,m_g,sigma_g,omega,ord); % self impedances of the overhead conductors
        Zm_over=Z_wise_mut(h,d,e_g,m_g,sigma_g,omega,ord); % mutual impedances of the overhead conductors
        kxm=0;
        Zm_new=Z_new_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxm); % mutual impedances in the mixed configuration
        % Total matrices
        Zg_OverUnder=Zs_papad+Zm_papad+Zs_over+Zm_over+Zm_new;
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
        Zm_pol_ovund=Z_pol_mut_overunder(h,d,sigma_g,f,ord); % mutual impedances in the mixed configuration Pollaczek
        Zs_pol_und=Z_pol_slf_und(h,cab_ex,sigma_g,f,ord);
        Zm_pol_und=Z_pol_mut_und(h,d,sigma_g,f,ord);
        Zs_pol_ov=Z_carson_slf(h,cab_ex,sigma_g,omega,0*e_g,ord);
        Zm_pol_ov=Z_carson_mut(h,d,omega,sigma_g,0*e_g,ord);
        % Total matrices
        Zg_CarsonPol=Zs_pol_und+Zm_pol_und+Zs_pol_ov+Zm_pol_ov+Zm_pol_ovund;
        Z_CarsonPol=Zin+Zg_CarsonPol;
        Ztot_CarsonPol(:,:,k) = bundleReduction(ph_order,Z_CarsonPol);
        if k==siz
            out(o).VarName='Ztot_CarsonPol';
            out(o).Label='Carson-Pollaczek';
            out(o).Values=Ztot_CarsonPol;
            o=o+1;
        end
    end
    
    %%%% De Conti (underground)
    if useFormula('Conti') && is_underground
        if k==1; Ztot_Conti=zeros(Nph,Nph,siz); end % Prelocate matrix
        Zs_conti=Z_xue_closed_slf(abs(h),cab_ex,e_g,sigma_g,f,ord);
        Zm_conti=Z_xue_closed_mut(abs(h),d,e_g,sigma_g,f,ord);
        % Total matrices
        Zg_conti=Zs_conti+Zm_conti;
        Z_Conti=Zin+Zg_conti;
        Ztot_Conti(:,:,k) = bundleReduction(ph_order,Z_Conti);
        if k==siz
            out(o).VarName='Ztot_Conti';
            out(o).Label='De Conti (underground)';
            out(o).Values=Ztot_Conti;
            o=o+1;
        end
    end

    %%%% FEMM solver
    if useFormula('FEMM')

        foldername=fullfile(currPath,'femm_files_Z');
        if isfolder(foldername)
            rmdir(foldername, 's');
        end
        mkdir(foldername)

        basename=fullfile(foldername,...
                sprintf('femmZ_%s_f%d',jobid,k));

        if k==1; Ztot_FEMM=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self and Mutual Impedances, seems to include skin effects
        Z_FEMM = Z_femm_slf_mut(Geom,soilFD,k,f,ord,basename);
        Ztot_FEMM(:,:,k) = bundleReduction(ph_order,Z_FEMM);
        % Store outputs
        if k==siz
            out(o).VarName='Ztot_FEMM';
            out(o).Label='FEMM';
            out(o).Values=Ztot_FEMM;
            o=o+1;
        end
    end
end % of main loop

% Add decent labels to uniform soil results if a layered structure is
% present
if soilFD.num_layers > 1
    for i = 1:numel(out)
        % Get the current label
        label = out(i).Label;
        
        % Check and skip if layered formulation 
        if contains(label, 'layers', 'IgnoreCase', true) || contains(label, 'FEMM', 'IgnoreCase', true)
            continue;  % Skip to the next iteration, leave this label untouched
        end

        % Check if the label already has something in parentheses at the end
        if regexp(label, '\(.*\)$')  % Regular expression to find existing parentheses at the end
            % Append tag inside the existing parentheses
            out(i).Label = regexprep(label, '\((.*)\)$', ['($1, ' append_tag ')']);
        else
            % No parentheses at the end, so just add tag at the end
            out(i).Label = [label ' (' append_tag ')'];
        end
    end
end

%% TESTME
if useFormula('ImportCYZ')
    fname=opts.CYZfile;
    if isfile(fname)
        [w, CYZ, ~] = get_ZY_from_cyz(fname);
        if size(CYZ,1) ~= Nph; error('The supplied LineCable_Data file has a different number of phases than the current model. Comparing oranges to potatoes eh?'); end
        CYZ = interp_matrix(CYZ, w, omega_total);
        out(o).VarName='Ztot_CYZ';
        out(o).Label='LineCableData';
        out(o).Values=CYZ;
        o=o+1;
    else
        disp('No valid CYZ file found in the job directory. Skipping...');
    end
end

if useFormula('ImportMAT')
    fname=opts.MATfile;
    if isfile(fname)
        [ff, Zmat, ~, varlbl] = get_ZY_from_mat(fname);
        if size(Zmat,1) ~= Nph; error('The supplied MAT-file has a different number of phases than the current model. Comparing oranges to apples eh?'); end
        Zmat = interp_matrix(Zmat, ff, f_total);
        out(o).VarName='Ztot_MAT';
        out(o).Label=varlbl;
        out(o).Values=Zmat;
        o=o+1;
    else
        disp('No valid MAT-file file found in the job directory. Skipping...');
    end
end

end % of main function

