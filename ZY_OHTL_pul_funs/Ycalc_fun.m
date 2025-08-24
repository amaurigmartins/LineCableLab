function [out] = Ycalc_fun(f_total,ord,Nph,soilFD,h,d,Geom,currPath,jobid,opts)

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

% Compute the actual parameters
e0=8.854187817e-12;  % Farads/meters
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
   
    %%%% Potential Coeficients of Conductor insulation
    Pin_mat=Pins_mat_fun(ord,Geom);
    Pin=Pin_mat./(2*pi*e0);

    %%%% Influence of perfect earth
    Ps_pet_perf=P_pet_slf_perf(h,cab_ex,ord);
    Pm_pet_perf=P_pet_mut_perf(h,d,ord);


    %%%% Pettersson
    if useFormula('Pettersson') && is_overhead
        if k==1; Ytot_Pet=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        Ps_pet_imperf=P_pet_slf_imperf(h,e_g,m_g,sigma_g,omega,ord);
        Ps_pet=Ps_pet_perf+Ps_pet_imperf;
        % Mutual
        Pm_pet_imperf=P_pet_mut_imperf(h,d,e_g,m_g,sigma_g,omega,ord);
        Pm_pet=Pm_pet_perf+Pm_pet_imperf;
        % Total matrices
        Pg_pet=Ps_pet+Pm_pet;
        P_pet=Pin_mat+Pg_pet;
        Ptot_pet = bundleReduction(ph_order,P_pet);
        Ytot_Pet(:,:,k)=1i.*omega.*e0.*2.*pi.*inv(Ptot_pet);
        % Store outputs
        if k==siz
            out(o).VarName='Ytot_Pet';
            out(o).Label='Pettersson';
            out(o).Values=Ytot_Pet;
            o=o+1;
        end
    end

    %%%% Images
    if useFormula('Images') && is_overhead
        if k==1; Ytot_Imag=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        Ps_imag=Ps_pet_perf;
        % Mutual
        Pm_imag=Pm_pet_perf;
        % Total matrices
        Pg_imag=Ps_imag+Pm_imag;
        P_imag=Pin_mat+Pg_imag;
        Ptot_imag = bundleReduction(ph_order,P_imag);
        Ytot_Imag(:,:,k)=1i.*omega.*e0.*2.*pi.*inv(Ptot_imag);
        % Store outputs
        if k==siz
            out(o).VarName='Ytot_Imag';
            out(o).Label='Images';
            out(o).Values=Ytot_Imag;
            o=o+1;
        end
    end

    %%%% Wise
    if useFormula('Wise') && is_overhead
        if k==1; Ytot_Wise=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        Ps_wise_imperf=P_wise_slf_imperf(h,e_g,m_g,sigma_g,omega,ord);
        Ps_wise=Ps_pet_perf+Ps_wise_imperf;
        % Mutual
        Pm_wise_imperf=P_wise_mut_imperf(h,d,e_g,m_g,sigma_g,omega,ord);
        Pm_wise=Pm_pet_perf+Pm_wise_imperf;
        % Total matrices
        Pg_wise=Ps_wise+Pm_wise;
        P_wise=Pin_mat+Pg_wise;
        Ptot_wise = bundleReduction(ph_order,P_wise);
        Ytot_Wise(:,:,k)=1i.*omega.*e0.*2.*pi.*inv(Ptot_wise);
        % Store outputs
        if k==siz
            out(o).VarName='Ytot_Wise';
            out(o).Label='Wise';
            out(o).Values=Ytot_Wise;
            o=o+1;
        end
    end

    % %%%% Kikuchi
    % if useFormula('Kikuchi') && is_overhead
    %     if k==1; Ytot_Kik=zeros(Nph,Nph,siz); end % Prelocate matrix
    %     % Self
    %     kxa='k0';
    %     Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
    %     % Mutual
    %     Pm_kik=P_kik_mut(h,d,e_g,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
    %     % Total matrices
    %     Pg_kik=Ps_kik+Pm_kik;
    %     P_kik=Pin+Pg_kik;
    %     Ptot_kik = bundleReduction(ph_order,P_kik);
    %     Ytot_Kik(:,:,k)=1i.*omega.*inv(Ptot_kik);
    %     % Store outputs
    %     if k==siz
    %         out(o).VarName='Ytot_Kik';
    %         out(o).Label='Kikuchi';
    %         out(o).Values=Ytot_Kik;
    %         o=o+1;
    %     end
    % end
    
    %%%% Papadopoulos (2-layered soil)
    if useFormula('Papad2Layers') && islayered_earth && is_overhead
        if k==1; Ytot_Papad2La=zeros(Nph,Nph,siz); end % Prelocate matrix
        % sigma_g_la=soilFD.sigma_g_total(1,:);
        % e_g_la=soilFD.e_g_total(1,:);
        % m_g_la=[soilFD.layer(:).m_g soilFD.m_g];
        kxa='k0';
        t=-soilFD.layer(1).t;
        % Self
        Ps2la=P_ohl_slf_2lay(h,cab_ex,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxa);
        % Mutual
        Pm2la=P_ohl_mut_2lay(h,d,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxa);
        % Total matrices
        Pg_Papad2La=Ps2la+Pm2la;
        P_Papad2La=Pin+Pg_Papad2La;
        Ptot_Papad2La = bundleReduction(ph_order,P_Papad2La);
        Ytot_Papad2La(:,:,k)=1i.*omega.*inv(Ptot_Papad2La);
        if k==siz
            out(o).VarName='Ytot_Papad2La';
            out(o).Label='Papadopoulos (overhead, 2-layers)';
            out(o).Values=Ytot_Papad2La;
            o=o+1;
        end
    end

    %%%% Papadopoulos (underground, 2-layered soil)
    if useFormula('Papad2LayersUnder') && islayered_earth && is_underground
        if k==1; Ytot_Papad2LaUnder=zeros(Nph,Nph,siz); end % Prelocate matrix
        % sigma_g_la=soilFD.sigma_g_total(1,:);
        % e_g_la=soilFD.e_g_total(1,:);
        % m_g_la=[soilFD.layer(:).m_g soilFD.m_g];
        kxe=1;
        t=-soilFD.layer(1).t;
        % Self
        Ps2la_under=P_papad_slf_2lay_under(h,cab_ex,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxe);
        % Mutual
        Pm2la_under=P_papad_mut_2lay_under(h,d,e_g_la,m_g_la,sigma_g_la,t,f,ord,kxe);
        % Total matrices
        Pg_Papad2LaUnder=Ps2la_under+Pm2la_under;
        P_Papad2LaUnder=Pin+Pg_Papad2LaUnder;
        Ptot_Papad2LaUnder = bundleReduction(ph_order,P_Papad2LaUnder);
        Ytot_Papad2LaUnder(:,:,k)=1i.*omega.*inv(Ptot_Papad2LaUnder);
        if k==siz
            out(o).VarName='Ytot_Papad2LaUnder';
            out(o).Label='Papadopoulos (underground, 2-layers)';
            out(o).Values=Ytot_Papad2LaUnder;
            o=o+1;
        end
    end
    
    %%%% Papadopoulos (underground)
    if useFormula('Papadopoulos') && is_underground
        if k==1; Ytot_Papad=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        kxe='k1';
        Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe);
        % Mutual
        Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe);
        % Total matrices
        Pg_Papad=Ps_papad+Pm_papad;
        P_Papad=Pin+Pg_Papad;
        Ptot_Papad = bundleReduction(ph_order,P_Papad);
        Ytot_Papad(:,:,k)=1i.*omega.*inv(Ptot_Papad);
        if k==siz
            out(o).VarName='Ytot_Papad';
            out(o).Label='Papadopoulos (underground)';
            out(o).Values=Ytot_Papad;
            o=o+1;
        end
    end

    %%%% Xue (underground)
    if useFormula('Xue') && is_underground
        if k==1; Ytot_Xue=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        Ps_xue=P_xue_slf(h,cab_ex,e_g,sigma_g,f,ord);
        % Mutual
        Pm_xue=P_xue_mut(h,d,e_g,sigma_g,f,ord);
        % Total matrices
        Pg_xue=Ps_xue+Pm_xue;
        P_Xue=Pin+Pg_xue;
        Ptot_Xue = bundleReduction(ph_order,P_Xue);
        Ytot_Xue(:,:,k)=1i.*omega.*inv(Ptot_Xue);
        if k==siz
            out(o).VarName='Ytot_Xue';
            out(o).Label='Xue (underground)';
            out(o).Values=Ytot_Xue;
            o=o+1;
        end
    end

    %%%% Martins-Papadopoulos-Chrysochos (overhead-underground)
    if useFormula('OverUnder')
        if k==1; Ytot_OverUnder=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Selfs and mutuals for underground
        kxe=0; %set 'k0' (string) to air layer, 'k1' for earth, 0 to neglect
        Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self coefficients of the underground conductors
        Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe); % mutual coefficients of the underground conductors
        % Selfs and mutuals for overhead
        kxa=0;
        % Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
        % Pm_kik=P_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        % TEMPORARY FALLBACK TO WISE
        Ps_over_imperf=P_wise_slf_imperf(h,e_g,m_g,sigma_g,omega,ord);
        Ps_over=(Ps_pet_perf+Ps_over_imperf)./(e0.*2.*pi);
        % Mutual
        Pm_over_imperf=P_wise_mut_imperf(h,d,e_g,m_g,sigma_g,omega,ord);
        Pm_over=(Pm_pet_perf+Pm_over_imperf)./(e0.*2.*pi);
        
        % Mutuals for overhead-underground
        kxm=0;
        Pm_overunder=P_new_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxm); % mutual coefficients in the mixed configuration
        % Total matrices
        Pg_OverUnder=Ps_papad+Pm_papad+Ps_over+Pm_over+Pm_overunder;
        P_OverUnder=Pin+Pg_OverUnder;
        Ptot_OverUnder = bundleReduction(ph_order,P_OverUnder);
        Ytot_OverUnder(:,:,k)=1i.*omega.*inv(Ptot_OverUnder);
        if k==siz
            out(o).VarName='Ytot_OverUnder';
            out(o).Label='Martins-Papadopoulos-Chrysochos (overhead-underground)';
            out(o).Values=Ytot_OverUnder;
            o=o+1;
        end
    end

    %%%% Martins-Papadopoulos-Chrysochos (no capacitive coupling)
    if useFormula('NoCapCoupling')
        if k==1; Ytot_NoCapCoupling=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Selfs and mutuals for underground
        kxe=0; %set 'k0' (string) to air layer, 'k1' for earth, 0 to neglect
        Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self coefficients of the underground conductors
        Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe); % mutual coefficients of the underground conductors
        % Selfs and mutuals for overhead
        kxa=0;
        % Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
        % Pm_kik=P_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        
        % TEMPORARY FALLBACK TO WISE
        Ps_over_imperf=P_wise_slf_imperf(h,e_g,m_g,sigma_g,omega,ord);
        Ps_over=(Ps_pet_perf+Ps_over_imperf)./(e0.*2.*pi);
        % Mutual
        Pm_over_imperf=P_wise_mut_imperf(h,d,e_g,m_g,sigma_g,omega,ord);
        Pm_over=(Pm_pet_perf+Pm_over_imperf)./(e0.*2.*pi);
        
        % Mutuals for overhead-underground
        kxm=0;
        Pm_overunder=zeros(ord,ord); % mutual coefficients in the mixed configuration
        % Total matrices
        Pg_NoCapCoupling=Ps_papad+Pm_papad+Ps_over+Pm_over+Pm_overunder;
        P_NoCapCoupling=Pin+Pg_NoCapCoupling;
        Ptot_NoCapCoupling = bundleReduction(ph_order,P_NoCapCoupling);
        Ytot_NoCapCoupling(:,:,k)=1i.*omega.*inv(Ptot_NoCapCoupling);
        if k==siz
            out(o).VarName='Ytot_NoCapCoupling';
            out(o).Label='Martins-Papadopoulos-Chrysochos (without C coupling)';
            out(o).Values=Ytot_NoCapCoupling;
            o=o+1;
        end
    end

    %%%% Overhead-underground arrangement, classical TL
    if useFormula('OverUnderClassicalTL')
        % if k==1; Ytot_OverUnderClassicalTL=zeros(Nph,Nph,siz); end % Prelocate matrix
        % if ~useFormula('OverUnder') % skip if already calculated
        %     % Selfs and mutuals for underground
        %     kxe=0; %set 'k0' (string) to air layer, 'k1' for earth, 0 to neglect
        %     Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self coefficients of the underground conductors
        %     Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe); % mutual coefficients of the underground conductors
        %     % Selfs and mutuals for overhead
        %     kxa=0;
        %     Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
        %     Pm_kik=P_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        % end
        % % Mutuals for overhead-underground are neglected
        % Pm_overunderClassicalTL=zeros(ord,ord);
        % % Total matrices
        % Pg_OverUnderClassicalTL=Ps_papad+Pm_papad+Ps_kik+Pm_kik+Pm_overunderClassicalTL;
        % P_OverUnderClassicalTL=Pin+Pg_OverUnderClassicalTL;
        % Ptot_OverUnderClassicalTL = bundleReduction(ph_order,P_OverUnderClassicalTL);
        % Ytot_OverUnderClassicalTL(:,:,k)=1i.*omega.*inv(Ptot_OverUnderClassicalTL);
        % if k==siz
        %     out(o).VarName='Ytot_OverUnderClassicalTL';
        %     out(o).Label='Martins-Papadopoulos-Chrysochos (earth shield)';
        %     out(o).Values=Ytot_OverUnderClassicalTL;
        %     o=o+1;
        % end
        if k==1; Ytot_OverUnderClassicalTL=zeros(Nph,Nph,siz); end % Prelocate matrix
        
        % Selfs and mutuals for underground
        Ps_papad=P_papad_slf(h,cab_ex,0*e_g,m_g,sigma_g,f,ord,0); % self coefficients of the underground conductors
        Pm_papad=P_papad_mut(h,d,     0*e_g,m_g,sigma_g,f,ord,0); % mutual coefficients of the underground conductors
        
        % Selfs and mutuals for overhead
        Ps_img=P_pet_slf_perf(h,cab_ex,ord)./(e0.*2.*pi); % True classical
        Pm_img=P_pet_mut_perf(h,d,ord)./(e0.*2.*pi);

        % Mutuals for overhead-underground are neglected
        Pm_overunderClassicalTL=zeros(ord,ord);
        
        % Total matrices
        Pg_OverUnderClassicalTL=Ps_papad+Pm_papad+Ps_img+Pm_img+Pm_overunderClassicalTL;
        P_OverUnderClassicalTL=Pin+Pg_OverUnderClassicalTL;
        Ptot_OverUnderClassicalTL = bundleReduction(ph_order,P_OverUnderClassicalTL);
        Ytot_OverUnderClassicalTL(:,:,k)=1i.*omega.*inv(Ptot_OverUnderClassicalTL);
        if k==siz
            out(o).VarName='Ytot_OverUnderClassicalTL';
            out(o).Label='Classical TL (overhead-underground)';
            out(o).Values=Ytot_OverUnderClassicalTL;
            o=o+1;
        end
    end

    %%%% De Conti (underground)
    if useFormula('Conti') && is_underground
        if k==1; Ytot_Conti=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        Ps_conti=P_xue_closed_slf(abs(h),cab_ex,e_g,sigma_g,f,ord);
        % Mutual
        Pm_conti=P_xue_closed_mut(abs(h),d,e_g,sigma_g,f,ord);
        % Total matrices
        Pg_conti=Ps_conti+Pm_conti;
        P_Conti=Pin+Pg_conti;
        Ptot_Conti = bundleReduction(ph_order,P_Conti);
        Ytot_Conti(:,:,k)=1i.*omega.*inv(Ptot_Conti);
        if k==siz
            out(o).VarName='Ytot_Conti';
            out(o).Label='De Conti (underground)';
            out(o).Values=Ytot_Conti;
            o=o+1;
        end
    end

    %%%% FEMM solver
    if useFormula('FEMM')

        foldername=fullfile(currPath,'femm_files_Y');
        if isfolder(foldername)
            rmdir(foldername, 's');
        end
        mkdir(foldername)

        basename=fullfile(foldername,...
                sprintf('femmY_%s_f%d',jobid,k));

        if k==1; Ytot_FEMM=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self and Mutual Admittances, seems to account for insulation
        Y_FEMM = Y_femm_slf_mut(Geom,soilFD,k,f,ord,basename);
        Ytot_FEMM(:,:,k) = bundleReduction(ph_order,Y_FEMM);
        % Store outputs
        if k==siz
            out(o).VarName='Ytot_FEMM';
            out(o).Label='FEMM';
            out(o).Values=Ytot_FEMM;
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
        [w,~ , CYZ] = get_ZY_from_cyz(fname);
        if size(CYZ,1) ~= Nph; error('The supplied LineCable_Data file has a different number of phases than the current model. Comparing oranges to potatoes eh?'); end
        CYZ = interp_matrix(CYZ, w, omega_total);
        out(o).VarName='Ytot_CYZ';
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
        [ff, ~, Ymat, varlbl] = get_ZY_from_mat(fname);
        if size(Ymat,1) ~= Nph; error('The supplied MAT-file has a different number of phases than the current model. Brilliant.'); end
        Ymat = interp_matrix(Ymat, ff, f_total);
        out(o).VarName='Ytot_MAT';
        out(o).Label=varlbl;
        out(o).Values=Ymat;
        o=o+1;
    else
        disp('No valid MAT-file file found in the job directory. Skipping...');
    end
end

end % of main function

