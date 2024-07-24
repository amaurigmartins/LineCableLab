function [] = makeJMartiModel_real(f,Z,Y,Ti_dis,g_dis,line_length,jobid,currPath,optsRW)

if ~isfield(optsRW,'plot1');optsRW.plot1=false;end
if ~isfield(optsRW,'plot2');optsRW.plot2=false;end
if ~isfield(optsRW,'Niter');optsRW.Niter = [20 50];end
if ~isfield(optsRW,'TendsToZero');optsRW.TendsToZero = true;end
if ~isfield(optsRW,'Npoles');optsRW.Npoles = [2 20];end
if ~isfield(optsRW,'samePolesH');optsRW.samePolesH = false;end

if size(f,1) == 1
    f=f.';
end

freq_siz=length(f);

%% Calculating the frequency-dependent  real transformation matrix

% Raw impedance and admittance matrices
if size(Z,3)==freq_siz
    Z=permute(Z,[3 1 2]);
    Y=permute(Y,[3 1 2]);
end

% Frequency-variant transform matrix from OHLT
ord=size(Z,2);
Ti=list2sqmat(Ti_dis,ord,freq_siz);

% Now let's pick the corresponding matrix for a specific frequency @ f_Ti
f_Ti=1e6;
interplM2freq = @(M) squeeze(interp1(f,M,f_Ti)); %use interpolation to handle the case where the frequency sample is missing
T=interplM2freq(Ti);
T=real(T); 

%% Calculating the characteristic impedance, the propagation function and phase velocity

% Recompute Zc and H using modified T and g
modif_T=sqmat2list(T,ord,freq_siz);
[Zch_m,~,~,~]=calc_char_imped_admit(modif_T,Z,Y,ord,freq_siz);
% Zch_m_OHLT=Zch_mod;
[H_m,~,~,~,~]=calc_prop_function(modif_T,g_dis,line_length,ord,freq_siz,f); %it works!!!

% Calculate phase velocities
for m=1:ord
    vel(:,m)=(2.*pi.*f)./imag(g_dis(:,m));
end


%% Calculating Zc poles and residues

if optsRW.samePolesH
        
    [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(Zch_m,f,optsRW); %Calculating the poles and residues
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
            [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(Zch_m(:,m),f,optsRW); %Calculating the poles and residues
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
        semilogx(f,abs(Zch_m(:,m)),'o', 'DisplayName', ['Mode #' num2str(m)]);hold all;
        semilogx(f,abs(fitOHLT_Zc(m).ffit), 'DisplayName', ['Fit mode #' num2str(m)]);hold all;
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
        [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(H_shifted,f,optsRW); %Calculating the poles
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
                [fit_data_JM, ffit_JM, numpol_JM,rmserr_JM ] = rationalfit_wrapper_real(H_shifted(:,m),f,optsRW); %Calculating the poles and residues
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
        semilogx(f,abs((H_m(:,m))),'o', 'DisplayName', ['Mode #' num2str(m)]);hold all;
        semilogx(f,abs(fitOHLT_H(m).ffit(1:freq_siz)), 'DisplayName', ['Fit mode #' num2str(m)]);hold all;
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
    tmin = 0.8*tmin;
end

options = optimset('Display','none', 'TolX',1e-6);
tau_opt=fminbnd(@(x)fcalc_WRP(H,f,x,optsRW),tmin,tmax,options); %Find optimal time delay
H_shifted = H.*exp(1i.*2.*pi.*f.*tau_opt);

end


function [rmserr_JM] = fcalc_WRP(H, f, tau,optsRW)
H_mod = H.*exp(1i.*2.*pi.*f.*tau);
[~, ~, ~,rmserr_JM ] = rationalfit_wrapper_real(H_mod,f,optsRW); %Find the RMSerr
end