function [] = makeJMartiModel_real(f,Z,Y,Ti_dis,g_dis,line_length,jobid,currPath,optsRW)

if ~isfield(optsRW,'plot1');optsRW.plot1=false;end
if ~isfield(optsRW,'plot2');optsRW.plot2=false;end
if ~isfield(optsRW,'Niter');optsRW.Niter = [20 50];end
if ~isfield(optsRW,'TendsToZero');optsRW.TendsToZero = true;end
if ~isfield(optsRW,'Npoles');optsRW.Npoles = [1 20];end
if ~isfield(optsRW,'samePolesH');optsRW.samePolesH = false;end
if ~isfield(optsRW,'tol');optsRW.tol = 3;end
if ~isfield(optsRW,'bode');optsRW.bode = false;end
if ~isfield(optsRW,'errtgt');optsRW.errtgt = 1e-3;end
if ~isfield(optsRW,'err_par');optsRW.err_par = true;end


if size(f,1) == 1
    f=f.';
end

freq_siz=length(f);

% Raw impedance and admittance matrices
if size(Z,3)==freq_siz
    Z=permute(Z,[3 1 2]);
    Y=permute(Y,[3 1 2]);
end

ord=size(Z,2);

%% Run optimization to find the best frequency sample for T

[H_test,~,~,~,~]=calc_prop_function(Ti_dis,g_dis,line_length,ord,freq_siz,f); 

if optsRW.samePolesH
            if optsRW.bode
                [~, ffit_JM,~,~] = rationalfit_wrapper_real_bode(H_test,f,optsRW); %Calculating the poles and residues
            else
                [~, ffit_JM,~,~] = rationalfit_wrapper_real(H_test,f,optsRW); %Calculating the poles and residues
            end
            for m=1:ord
                H_modefit(:,m) = ffit_JM(:,m);
            end

else
            for m=1:ord
                if optsRW.bode
                    [~, ffit_JM,~,~] = rationalfit_wrapper_real_bode(H_test(:,m),f,optsRW); %Calculating the poles and residues
                else
                    [~, ffit_JM,~,~] = rationalfit_wrapper_real(H_test(:,m),f,optsRW);
                end
                H_modefit(:,m) = ffit_JM;
            end
            
end
T_const = real(Ti_dis);
T_const=list2sqmat(T_const,ord,freq_siz);

for k=1:freq_siz
    for y=1:freq_siz
        X(:,:,y) = H_test(y,:) - ((squeeze(T_const(k,:,:)))*(diag(H_modefit(y,:)))*(squeeze(T_const(k,:,:)))); 
        err_T_const(y,1,k) = norm(X(:,:,y),"fro");
    end
    [err_T(k,1),p_T(k,1)] = min(err_T_const(:,1,k));
end

% Now let's pick the corresponding matrix for a specific frequency @ f_Ti
[~,p_Ti] = min(err_T);
f_Ti=f(p_T(p_Ti,1));
interplM2freq = @(M) squeeze(interp1(f,M,f_Ti)); %use interpolation to handle the case where the frequency sample is missing
T=interplM2freq(T_const);
% T=real(T); 

%% Calculating the characteristic impedance, the propagation function and phase velocity

% Zc and H using the chosen T 
modif_T_Jm=sqmat2list(T,ord,freq_siz);
[Zch_m,~,~,~]=calc_char_imped_admit(modif_T_Jm,Z,Y,ord,freq_siz);
[H_m,~,~,~,~]=calc_prop_function(modif_T_Jm,g_dis,line_length,ord,freq_siz,f); %it works!!!

% Calculate phase velocities
for m=1:ord
    vel(:,m)=(2.*pi.*f)./imag(g_dis(:,m));
end


%% Calculating Zc poles and residues
if optsRW.samePolesH
    if optsRW.bode   
        [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real_bode(Zch_m,f,optsRW); %Calculating the poles and residues
    else
        [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(Zch_m,f,optsRW); %Calculating the poles and residues
    end
        for m=1:ord
            fitOHLT_Zc(m).mode = m;
            fitOHLT_Zc(m).NORD = numpol_JM;
            fitOHLT_Zc(m).ks = fit_data_JM(m).D;
            fitOHLT_Zc(m).pol = -fit_data_JM(m).A;
            fitOHLT_Zc(m).res = fit_data_JM(m).C;
            fitOHLT_Zc(m).err = rmserr_JM;
            fitOHLT_Zc(m).ffit = ffit_JM(:,m);
        end
else
        for m=1:ord
            if optsRW.bode
                [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real_bode(Zch_m(:,m),f,optsRW); %Calculating the poles and residues
            else
                [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(Zch_m(:,m),f,optsRW); %Calculating the poles and residues
            end
            fitOHLT_Zc(m).mode = m;
            fitOHLT_Zc(m).NORD = numpol_JM;
            fitOHLT_Zc(m).ks = fit_data_JM.D;
            fitOHLT_Zc(m).pol = -fit_data_JM.A;
            fitOHLT_Zc(m).res = fit_data_JM.C;
            fitOHLT_Zc(m).err = rmserr_JM;
            fitOHLT_Zc(m).ffit = ffit_JM;
        end
     
end


% Plot for Zc fitting x data (N_phases)
if optsRW.plot1
    figure
    for m=1:ord
        semilogx(f,abs(Zch_m(:,m)), 'DisplayName', ['Mode #' num2str(m)]);hold all;
        semilogx(f,abs(fitOHLT_Zc(m).ffit),'o', 'DisplayName', ['Fit mode #' num2str(m)]);hold all;
    end
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('Z_{ch} [\Omega]')
    grid on;grid minor;
    legend('Location','best')
end

%% Calculating H with optimal time delay
for m=1:ord
    [H_shifted_,tau_opt]=findoptimtau_WRP(f,vel(:,m),H_m(:,m),line_length,optsRW); % Optimized time finding function
    Optimal_time(m,1) = tau_opt;
    H_shifted(:,m)= H_shifted_;
end

%% Calculating fitting of H for the modes
if optsRW.samePolesH
    if optsRW.bode
        [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real_bode(H_shifted,f,optsRW); %Calculating the poles
    else
        [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(H_shifted,f,optsRW); %Calculating the poles
    end
            for m=1:ord
                fitOHLT_H(m).mode = m;
                fitOHLT_H(m).NORD = numpol_JM;
                fitOHLT_H(m).tau_opt = Optimal_time(m,1);
                fitOHLT_H(m).pol = -fit_data_JM(m).A;
                fitOHLT_H(m).res = fit_data_JM(m).C;
                fitOHLT_H(m).err = rmserr_JM;
                fitOHLT_H(m).ffit = ffit_JM(:,m);
            end

else
            for m=1:ord
                if optsRW.bode
                    [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real_bode(H_shifted(:,m),f,optsRW); %Calculating the poles and residues
                else
                    [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(H_shifted(:,m),f,optsRW); %Calculating the poles and residues
                end
                fitOHLT_H(m).mode = m;
                fitOHLT_H(m).NORD = numpol_JM;
                fitOHLT_H(m).tau_opt = Optimal_time(m,1);
                fitOHLT_H(m).pol = -fit_data_JM.A;
                fitOHLT_H(m).res = fit_data_JM.C;
                fitOHLT_H(m).err = rmserr_JM;
                fitOHLT_H(m).ffit = ffit_JM;
            end
end
% Plot for H fitting x data (N_phases)

if optsRW.plot1
    figure
    for m=1:ord
        semilogx(f,abs((H_m(:,m))), 'DisplayName', ['Mode #' num2str(m)]);hold all;
        semilogx(f,abs(fitOHLT_H(m).ffit(1:freq_siz)),'o', 'DisplayName', ['Fit mode #' num2str(m)]);hold all;
    end
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('H')
    grid on
    legend ('Location','best')
end

% Plot for each mode of H and Zc

if optsRW.plot1
    for m=1:ord
        figure
        semilogx(f,abs(Zch_m(:,m)),'LineWidth',2);hold all; semilogx(f,abs(fitOHLT_Zc(m).ffit),'-.','LineWidth',2)
        legend('Data', 'Fit')
        title(sprintf('Zc mode #%d',m))
        axis tight
        xlabel('Frequency [Hz]')
        ylabel('Z_{ch} [\Omega]')
        grid on
        legend('Location','best')

        figure;semilogx(f,abs(H_m(:,m)),'LineWidth',2);hold all; semilogx(f,abs(fitOHLT_H(m).ffit),'-.','LineWidth',2)
        legend('Data', 'Fit')
        title(sprintf('H mode #%d',m))
        axis tight
        xlabel('Frequency [Hz]')
        ylabel('H')
        grid on
        legend('Location','best')
    end
end

%% Write data to PCH file
pchfname = fullfile(currPath,['PCH_JM_' jobid '.pch']);
fcontent = punchJMartiCard(ord, fitOHLT_Zc, fitOHLT_H, T, pchfname);

end

function [H_shifted,tau_opt] = findoptimtau_WRP(f,vel,H,line_length,optsRW)
j=length(f)-1;
w=2*pi*f;
Ns=length(w);
absH=abs(H);

phase1=(pi/2)*log((absH(j+1)/absH(j-1)))/(log(w(j+1)/w(j-1)));
phase2=0;
term2=log((absH(j+1)/absH(j-1))) /(log(w(j+1)/w(j-1)));
for k=2:Ns-1
    term1=log(absH(k+1)/absH(k-1)) /(log(w(k+1)/w(k-1)));
    if k~=j
        phase2=phase2+(abs(term1)-abs(term2))*log(coth(abs(log(w(k)/w(j)))/2))*log(w(k+1)/w(k));
    end
end

phase2=phase2/pi;
ph=(phase1-phase2);  
tau_mps=(line_length/vel(j)) + ph/(w(j));

tmin=1*abs(tau_mps);
tmax=1*abs(line_length/vel(end));

if tmin >= tmax
    tau_a=line_length/3e8;
    tau_b=(line_length/vel(j));
    tmin=.8*tau_a;
    tmax=1.2*tau_b;
end

options = optimset('Display','none', 'TolX',1e-6);
tau_opt=fminbnd(@(x)fcalc_WRP(H,f,x,optsRW),tmin,tmax,options); %Find optimal time delay
H_shifted = H.*exp(1i.*2.*pi.*f.*tau_opt);

end


function [rmserr_JM] = fcalc_WRP(H, f, tau,optsRW)
    H_mod = H.*exp(1i.*2.*pi.*f.*tau);
    if optsRW.bode
        [~, ~, ~,rmserr_JM ] = rationalfit_wrapper_real_bode(H_mod,f,optsRW); %Find the RMSerr
    else
        [~, ~, ~,rmserr_JM ] = rationalfit_wrapper_real(H_mod,f,optsRW); %Find the RMSerr
    end
end
