
ZYprnt=false;
useBodeFit=true;

% General parameters
f_Ti=1e6; % Frequency chosen for JMarti model
Ns = length(f); % Number of samples
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency
f = f.';

% Raw impedance and admittance matrices
Z=permute(Z,[3 1 2]);
Y=permute(Y,[3 1 2]);

% Frequency-variant transform matrix from OHLT
ord=size(Z,2);
Ti=list2sqmat(Ti_dis,ord,freq_siz);

% Now let's pick the corresponding matrix for a specific frequency @ f_Ti
interplM2freq = @(M) squeeze(interp1(f,M,f_Ti)); %use interpolation to handle the case where the frequency sample is missing
T=interplM2freq(Ti);
T=real(T); %is this really necessary?

% Recompute Zc and H using modified T and g
modif_T=sqmat2list(T,ord,freq_siz);
[Zch_m,Ych_m,~,~]=calc_char_imped_admit(modif_T,Z,Y,ord,freq_siz);
Zch_m_OHLT=Zch_mod;
[H_m,~,~,~,~]=calc_prop_function(modif_T,g_dis,line_length,ord,freq_siz,f); %it works!!!

% Calculate phase velocities
for m=1:ord
    vel(:,m)=(2.*pi.*f)./imag(g_dis(:,m));
end

% standard VF options
opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles 
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %Do NOT skip identification of residues (C,D,E) 
opts.cmplx_ss=1;   %Create complex state space model
opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=0;       %Create magnitude plot for fitting of f(s) 
opts.logx=1;       %Use logarithmic abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis 
opts.errplot=1;    %Include deviation in magnitude plot
opts.phaseplot=0;  %Also produce plot of phase angle (in addition to magnitiude)
opts.legend=1;     %Do include legends in plots

% options used by VF wrapper
opts.firstguesstype=1; %real valued poles
opts.output_messages=false;
opts.passive=0;

% VF impedance characteristic 
opts.weightscheme=1;  %1=uniform weights; 2=inv(sqrt(norm)) 3=inv(norm)
opts.asymp=2;      %%1-> d=e=0, 2 -> d!=0, e=0, 3-> d!=0, e!=0
Zch_src=Zch_m;
ERR=.1/100;
for m=1:ord
    fun=(Zch_src(:,m));
    [pol, res, ks, NORD, ffit, err] = vectfit_wrapper(fun,f,ERR,opts);
    if ~isreal(pol)
        warning(sprintf('VF was unable to fit the characteristic impedance Zc for mode #%d with real poles only. ATP results might be unstable or inaccurate.',m));
        if useBodeFit
            warning(sprintf('Will try to fit Zc for mode #%d using Bode method. Beware possible innaccuracies.',m));
            tol = 3; % Tolerance in decibels to set a new pole and/or new zero
            [P,Z,k] = Bode_process(abs(fun),f,Ns,tol); % Bode subroutine
            as = k.*poly(Z); bs = poly(P); % Polynomials
            [r,p,ks] = residue(as,bs); % Poles, residues and constant term
            TF=isempty(ks);if(TF==1);ks=0;end

            ffit = zeros(1,Ns)';
            for k = 1:length(p)
                ffit = ffit + (r(k)./(s.' - p(k))).';
            end
            ffit = ffit + ks;
            error = abs(fun - ffit);

            NORD = length(p);
            pol = -p;
            res = r;
        end
    end
    fitOHLT_Zc(m).mode = m;
    fitOHLT_Zc(m).NORD = NORD;
    fitOHLT_Zc(m).ks = ks;
    fitOHLT_Zc(m).pol = -pol;
    fitOHLT_Zc(m).res = res;
    fitOHLT_Zc(m).err = err;
    fitOHLT_Zc(m).ffit = ffit;
end

if ZYprnt
    %%% DEBUGME
    figure
    for m=1:ord
        semilogx(f,abs(Zch_src(:,m)),'o', 'DisplayName', ['mode #' num2str(m)]);hold all;
        semilogx(f,abs(fitOHLT_Zc(m).ffit), 'DisplayName', ['fit mode #' num2str(m)]);hold all;
    end
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('Z_{ch} [\Omega]')
    grid on
    legend
end

% VF propagation function with optimum delay time
opts.asymp=1;      %%1-> d=e=0, 2 -> d!=0, e=0, 3-> d!=0, e!=0
opts.weightscheme=3; %1=uniform weights; 2=inv(sqrt(norm)) 3=inv(norm)
H_src=H_m;
ERR=.1/100;
for m=1:ord
    [err, pol, res, ks, NORD, ffit, tau_mps, tau_opt, tau]=findoptimtau(f,vel(:,m),H_src(:,m),line_length,ERR,opts);

    if ~isreal(pol)
        warning(sprintf('VF was unable to fit the propagation function H for mode #%d with real poles only. Results in ATP may be unstable and/or overflow.',m));
        if useBodeFit
            warning(sprintf('Will try to fit A1 for mode #%d using Bode method. Beware possible innaccuracies.',m));
            A1min=.1;
            tol = 3; % Tolerance in decibels to set a new pole and/or new zero
            [ff, H] = fitnextrap(f,H,ord, 'decdecay');
            ss=1j.*2.*pi.*ff;
            A1 = H(:,m).*exp(1i*2*pi.*ff.*tau_opt(m));
            fun=abs(A1);
            [P,Z,k] = Bode_process(fun(fun>=A1min),ff(fun>=A1min),length(ff(fun>=A1min)),tol); % Bode subroutine
            as = k.*poly(Z); bs = poly(P); % Polynomials
            [r,p,ks] = residue(as,bs); % Poles, residues and constant term
            TF=isempty(ks);if(TF==1);ks=0;end
            ffit = zeros(1,length(ff))';
            for k = 1:length(p)
                ffit = ffit + (r(k)./(ss.' - p(k))).';
            end
            ffit = ffit + ks;

            error = abs(fun(1:Ns) - ffit(1:Ns));

            NORD = length(p);
            pol = -p;
            res = r;
        end
    end
    fitOHLT_H(m).mode = m;
    fitOHLT_H(m).NORD = NORD;
    fitOHLT_H(m).tau_opt = tau;
    fitOHLT_H(m).pol = -pol;
    fitOHLT_H(m).res = res;
    fitOHLT_H(m).err = err;
    fitOHLT_H(m).ffit = ffit;
end

if ZYprnt
    %%% DEBUGME
    figure
    for m=1:ord
        semilogx(f,abs((H_src(:,m))),'o', 'DisplayName', ['mode #' num2str(m)]);hold all;
        semilogx(f,abs(fitOHLT_H(m).ffit(1:freq_siz)), 'DisplayName', ['fit mode #' num2str(m)]);hold all;
    end
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('H')
    grid on
    legend
end

if ZYprnt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%DEBUGME
    for m=1:ord
        figure
        semilogx(f,abs(Zch_src(:,m)));hold all; semilogx(f,abs(Zch_m_OHLT(:,m)),'-.')
        legend('Const+real T @ chosen freq', 'Var Ti')
        title(sprintf('Zc mode #%d',m))
        axis tight
        xlabel('Frequency [Hz]')
        ylabel('Z_{ch} [\Omega]')
        grid on
        legend

        figure;semilogx(f,abs(H_src(:,m)));hold all; semilogx(f,abs(H_m_OHLT(:,m)),'-.')
        legend('Const+real T @ chosen freq', 'Var Ti')
        title(sprintf('H mode #%d',m))
        axis tight
        xlabel('Frequency [Hz]')
        ylabel('H')
        grid on
        legend
    end
end

% Now write data to PCH file
pchfname = fullfile(currPath,['pch_export__' jobid '.pch']);
fcontent = punchJMartiCard(ord, fitOHLT_Zc, fitOHLT_H, T, pchfname);