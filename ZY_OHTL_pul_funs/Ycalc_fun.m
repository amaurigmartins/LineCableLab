function [out] = Ycalc_fun(f_total,ord,Nph,soilFD,h,d,Geom,jobid,currPath,opts)

if nargin == 9; opts=struct();end

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

% Compute the actual parameters
e0=8.854187817e-12;  % Farads/meters
o=1; % counter for the number of outputs
for k=1:siz
    omega=omega_total(k);
    f=f_total(k);
    sigma_g=soilFD.sigma_g_total(k);
    m_g=soilFD.m_g(end);
    e_g=soilFD.e_g_total(k);
    
    %%%% Potential Coeficients of Conductor insulation
    Pin_mat=Pins_mat_fun(ord,Geom);
    Pin=Pin_mat./(2*pi*e0);

    %%%% Influence of perfect earth
    Ps_pet_perf=P_pet_slf_perf(h,cab_ex,ord);
    Pm_pet_perf=P_pet_mut_perf(h,d,ord);


    %%%% Pettersson
    if useFormula('Pettersson')
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
    if useFormula('Images')
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
    if useFormula('Wise')
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

    %%%% Kikuchi
    if useFormula('Kikuchi')
        if k==1; Ytot_Kik=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        global kxa;if isempty(kxa);kxa='k0';end;
        Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
        % Mutual
        Pm_kik=P_kik_mut(h,d,e_g,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        % Total matrices
        Pg_kik=Ps_kik+Pm_kik;
        P_kik=Pin+Pg_kik;
        Ptot_kik = bundleReduction(ph_order,P_kik);
        Ytot_Kik(:,:,k)=1i.*omega.*inv(Ptot_kik);
        % Store outputs
        if k==siz
            out(o).VarName='Ytot_Kik';
            out(o).Label='Kikuchi';
            out(o).Values=Ytot_Kik;
            o=o+1;
        end
    end
    
    %%%% Sunde (2-layered soil)
    if useFormula('Sunde2Layers')
        if k==1; Ytot_Sunde2La=zeros(Nph,Nph,siz); end % Prelocate matrix
        sigma_g_la=soilFD.sigma_g_total(1,:);
        e_g_la=soilFD.e_g_total(1,:);
        global kxa;if isempty(kxa);kxa='k0';end;
        t=-soilFD.layer(1).t;
        % Self
        Ps2la=P_ohl_slf_2lay(h,cab_ex,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        % Mutual
        Pm2la=P_ohl_mut_2lay(h,d,e_g_la,m_g,sigma_g_la,t,f,ord,kxa);
        % Total matrices
        Pg_Sunde2La=Ps2la+Pm2la;
        P_Sunde2La=Pin+Pg_Sunde2La;
        Ptot_Sunde2La = bundleReduction(ph_order,P_Sunde2La);
        Ytot_Sunde2La(:,:,k)=1i.*omega.*inv(Ptot_Sunde2La);
        if k==siz
            out(o).VarName='Ytot_Sunde2La';
            out(o).Label='Sunde (2-layers)';
            out(o).Values=Ytot_Sunde2La;
            o=o+1;
        end
    end
    
    %%%% Papadopoulos (underground)
    if useFormula('Papadopoulos')
        if k==1; Ytot_Papad=zeros(Nph,Nph,siz); end % Prelocate matrix
        % Self
        global kxe;if isempty(kxe);kxe='k1';end;
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
    if useFormula('Xue')
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
        global kxe;if isempty(kxe);kxe=0;end; %set 'k0' (string) to air layer, 'k1' for earth, 0 to neglect
        Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self coefficients of the underground conductors
        Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe); % mutual coefficients of the underground conductors
        % Selfs and mutuals for overhead
        global kxa;if isempty(kxa);kxa=0;end;
        Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
        Pm_kik=P_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        % Mutuals for overhead-underground
        global kxm;if isempty(kxm);kxm=0;end;
        Pm_overunder=P_new_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxm); % mutual coefficients in the mixed configuration
        % Total matrices
        Pg_OverUnder=Ps_papad+Pm_papad+Ps_kik+Pm_kik+Pm_overunder;
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

    %%%% Martins-Papadopoulos-Chrysochos (overhead-underground, classical TL)
    if useFormula('OverUnderClassicalTL')
        if k==1; Ytot_OverUnderClassicalTL=zeros(Nph,Nph,siz); end % Prelocate matrix
        if ~useFormula('OverUnder') % skip if already calculated
            % Selfs and mutuals for underground
            global kxe;if isempty(kxe);kxe=0;end; %set 'k0' (string) to air layer, 'k1' for earth, 0 to neglect
            Ps_papad=P_papad_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxe); % self coefficients of the underground conductors
            Pm_papad=P_papad_mut(h,d,e_g,m_g,sigma_g,f,ord,kxe); % mutual coefficients of the underground conductors
            % Selfs and mutuals for overhead
            global kxa;if isempty(kxa);kxa=0;end;
            Ps_kik=P_kik_slf(h,cab_ex,e_g ,m_g,sigma_g,f,ord,kxa); % self coefficients of the overhead conductors
            Pm_kik=P_kik_mut(h,d,e_g ,m_g,sigma_g,f,ord,kxa); % mutual coefficients of the overhead conductors
        end
        % Mutuals for overhead-underground are neglected
        Pm_overunderClassicalTL=zeros(ord,ord);
        % Total matrices
        Pg_OverUnderClassicalTL=Ps_papad+Pm_papad+Ps_kik+Pm_kik+Pm_overunderClassicalTL;
        P_OverUnderClassicalTL=Pin+Pg_OverUnderClassicalTL;
        Ptot_OverUnderClassicalTL = bundleReduction(ph_order,P_OverUnderClassicalTL);
        Ytot_OverUnderClassicalTL(:,:,k)=1i.*omega.*inv(Ptot_OverUnderClassicalTL);
        if k==siz
            out(o).VarName='Ytot_OverUnderClassicalTL';
            out(o).Label='Martins-Papadopoulos-Chrysochos (classical TL)';
            out(o).Values=Ytot_OverUnderClassicalTL;
            o=o+1;
        end
    end

end % of main loop

%% TESTME
if useFormula('CYZ')
    fname=fullfile(currPath,[jobid '.cyz']);
    if isfile(fname)
        [w,~ , CYZ] = get_ZY_from_cyz(fname);
        if size(CYZ,1) ~= Nph; error('The supplied LineCable_Data file has a different number of phases than the current model. Comparing oranges to potatoes eh?'); end
        CYZ = interp_matrix(CYZ, w, omega_total);
        out(o).VarName='Ytot_CYZ';
        out(o).Label='LineCable_Data';
        out(o).Values=CYZ;
        o=o+1;
    else
        disp('No valid CYZ file found in the job directory. Skipping...');
    end
end


end % of main function

