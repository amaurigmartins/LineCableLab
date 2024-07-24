function [] = makeULModel3(Z,Y,f,line_length,modif_T,g_dis,optsRW,jobid,currPath)


if ~isfield(optsRW,'plot1');optsRW.plot1=false;end;
if ~isfield(optsRW,'plot2');optsRW.plot2=false;end;
if ~isfield(optsRW,'Niter');optsRW.Niter = [20 50];end;
if ~isfield(optsRW,'TendsToZero');optsRW.TendsToZero = true;end;
if ~isfield(optsRW,'Npoles');optsRW.Npoles = [2 20];end;
if ~isfield(optsRW,'samePolesH');optsRW.samePolesH = true;end;

% Fit_print= false;    %Flag for plot for all fits
transp=@(x) x.';

%% Initial data

if size(f,1) == 1
    frequency=transp(f);
else
    frequency=f;
end

w = 2*pi.*frequency;    % Calculation of angular frequency (rad/sec)
s = transp(1j.*w);      % Calculation of the complex angular frequency
freq_siz = size(frequency,1);   % Number of frequency samples
ord=size(Z,2);  % Order of the system

%% Calculating the frequency-dependent transformation matrix

% Raw impedance and admittance matrices
if size(Z,3)==freq_siz
    Z=permute(Z,[3 1 2]);
    Y=permute(Y,[3 1 2]);
end

%  [modif_T,g_dis]=LM_calc_norm_str(ord,freq_siz,Z,Y,frequency); % Transformation matrix using the Levenberg-Marquardt (LM)

Ti = zeros(size(Z,3), size(Z,3), size(frequency,1));
invTi = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1:freq_siz
    for o = 1:ord
        Ti(o,:,k) = modif_T(k,(o-1)*ord+1:o*ord); % Ti in N_phases x N_phases x N_samples
    end
    invTi(:,:,k) = inv(Ti(:,:,k)); % Calculation the inverse of Ti
end

%% Calculating the characteristic impedance, admittance and the trace your trace

[~,~,~,Ych] = calc_char_imped_admit(modif_T,Z,Y,ord,freq_siz); % Calculation of characteristic impedance and admittance

Yc = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1:freq_siz
    for o = 1:ord
        Yc(o,:,k) = Ych(k,(o-1)*ord+1:o*ord); %  Ych in N_phases x N_phases x N_phases
    end
end

traceYc = zeros(size(frequency));
Yc_t = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1 : size(Z,1)
    Yc_t(:,:,k) = transp(Yc(:,:,k));
    traceYc(k) = trace( Yc(:,:,k) ); % Calculating the Yc trace
end

%% Calculating the propagation function and phase velocity
[Aj] = calc_prop_fun(modif_T,g_dis,line_length,ord,freq_siz); 	% Propagation function for each modes

vel = zeros(size(frequency,1), size(Z,3));
for m = 1:ord
    vel(:,m) = (2*pi*frequency)./imag(g_dis(:,m)); 	% Calculation of phase velocity
end
%% Calculating Yc trace poles
fitYc = transp(traceYc);
[fit_data_Yc, ffit, numpolYc, ~] = rationalfit_wrapper(fitYc,frequency,optsRW); %Calculating the poles

if ~ispassive(fit_data_Yc)
    warning('Resulting fitted traceYc is not passive. Results might be numerically unstable or inaccurate.')
end

if optsRW.plot1
    fig=1;
    figure(fig);
    semilogx(frequency,abs(traceYc),'-b', 'linewidth', 3); hold on
    semilogx(frequency,abs(ffit),'--r', 'linewidth', 3)
    legend('trace Yc', 'Fit trace Yc', 'location', 'north' )
    ylabel('Magnitude')
    xlabel('Frequency [Hz]')
    grid on
end
%% Calculating the residuals of Yc

% Opts for Vector Fitting
clear opts
opts.relax = 1;
opts.stable = 1;
opts.asymp = 2;
opts.skip_pole = 1;
opts.skip_res = 0;
opts.cmplx_ss = 1;
opts.spy1 = 0;
opts.spy2 = 0;

rYc = reshape(Yc_t, size(Yc_t,2)*size(Yc_t,1), size(Yc_t,3)); % Yc_t in  N_phases.N_phases x N_samples

for j = 1:1 % número de iterações
    weight = ones(1,numel(s)); % peso - 1 x Namostras
    [SERYc,~,~,Ycapprox,~] = vectfit3(rYc, s, transp(fit_data_Yc.A), weight, opts); % Function VF by Gustavsen
end

Ycapprox_ = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1:freq_siz
    for o = 1:ord
        Ycapprox_(o,:,k) = Ycapprox((o-1)*ord+1:o*ord,k); % Ycapprox in N_phases x N_phases x N_samples
    end
end

% Plot for each element of the Yc matrix (N_phases x N_phases)
if optsRW.plot1
    linewidth = 3;
    approxColor = ':r';
    for i = 1:ord
        for j = 1:ord
            fig_idx=fig+j;
            data = squeeze(Ych(:,i*j));
            approxData = squeeze(Ycapprox(i*j, :));
            figure(fig_idx); clf;
            for subp = 1:2
                if subp == 1
                    subplot(2, 1, subp)
                    semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Yc_{%d,%d}', i, j)); hold all;
                    semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpolYc)); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                else
                    subplot(2, 1, subp)
                    semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Yc_{%d,%d}', i, j)); hold all;
                    semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpolYc)); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
                end
            end
        end
        fig=fig_idx;
    end
end
%% Calculating Pj with optimal time delay

for m = 1:ord
    [tau_opt, Aj_shifted] = findoptimtau(frequency,vel(:,m),Aj(:,m),line_length,optsRW); 	% Optimized time finding function
    Optimal_time(m,1) = tau_opt;
    Pj(m,:)= Aj_shifted;
end

P_j_t = zeros(length(frequency),ord);
for i = 1:ord
    for j = 1:length(frequency)
        P_j_t(j,i) = transp(Pj(i,j));
    end
end

Pjj = zeros(size(P_j_t,2),size(P_j_t,2),size(f,2));
for k=1:freq_siz
    Pjj(:,:,k)=diag(P_j_t(k,:)); % Pj in a diagonal matrix
end

%% Calculating fitting of Pj with different poles for the modes
if optsRW.samePolesH % use the same set of poles to fit all modes of H
    fitPj_s=transp(Pj);
    [fit_data_s, ffit_s, ~, ~] = rationalfit_wrapper(fitPj_s,frequency,optsRW);	%Calculating the poles for each mode
    poles_wrapper = {fit_data_s.A};
    Pjapprox_s = ffit_s;
    %	Plot for the Pj fitting poles
    if optsRW.plot1
        fig=1;
        linewidth = 3;
        approxColor = ':r';
        for i = 1:ord
            fig_idx=i+fig;
            data = Pj(i,:);
            approxData = Pjapprox_s(:,i);
            figure(fig_idx); clf;
            for subp = 1:2
                if subp == 1
                    subplot(2, 1, subp)
                    semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Pj_{mode %d}', i)); hold all;
                    semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', size(poles_wrapper{1},1))); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                else
                    subplot(2, 1, subp)
                    semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Pj_{mode %d}', i)); hold all;
                    semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', size(poles_wrapper{1},1))); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
                end
            end
        end
        fig=fig_idx;
    end
else % each mode of H gets its own set of poles. Why? BECAUSE IT CAN BE DONE
    for t=1:size(Pj,1)
        fitPj=transp(Pj(t,:));
        [fit_data, ffit, ~, ~] = rationalfit_wrapper(fitPj,frequency,optsRW);	%Calculating the poles for each mode
        poles_wrapper(t) = {fit_data.A};
        Pjapprox(t,:) = ffit;
    end
    %	Plot for the Pj fitting poles
    if optsRW.plot1
        fig=1;
        linewidth = 3;
        approxColor = ':r';
        for i = 1:ord
            fig_idx=i+fig;
            data = Pj(i,:);
            approxData = Pjapprox(i,:);
            figure(fig_idx); clf;
            for subp = 1:2
                if subp == 1
                    subplot(2, 1, subp)
                    semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Pj_{mode %d}', i)); hold all;
                    semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', size(poles_wrapper{1},1))); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                else
                    subplot(2, 1, subp)
                    semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Pj_{mode %d}', i)); hold all;
                    semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', size(poles_wrapper{1},1))); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
                end
            end
        end
        fig=fig_idx;
    end
end

%% Calculating the residues of DjPj
for m=1:ord
    for k=1:freq_siz
        mode(m).D(:,:,k)=Ti(:,m,k)*invTi(m,:,k); % Idempotent Dj matrix -> Eliminates indeterminacy problems
    end
end

%	Opts for Vector fitting
clear opts
opts.relax = 1;
opts.stable = 1;
opts.asymp = 1;
opts.skip_pole = 1;
opts.skip_res = 0;
opts.cmplx_ss = 1;
opts.spy1 = 0;
opts.spy2 = 0;
opts.passive = 0;

for o = 1:ord
    clear DjPj
    for k = 1:freq_siz
        DjPj(:,:,k) = mode(o).D(:,:,k)*Pjj(o,o,k);
        DjPj_t(:,:,k) = transp(DjPj(:,:,k));
    end
    DjPj_ = reshape(DjPj_t, size(DjPj_t,2)*size(DjPj_t,1), size(DjPj_t,3)); %Yc_t in (N_phases*N_phases) x N_samples
    poles = poles_wrapper{o};
    numpoles_thismode = length(poles);
    for j = 1:1
        weight = ones(1,numel(s));
        [SER_DjPj,~,~,Happrox,~] = vectfit3(DjPj_, s, poles, weight, opts); %
        rDjPj(o) = {SER_DjPj.C};	% Residues for DjPj

        if optsRW.plot1
            linewidth = 3;
            approxColor = ':r';
            for i = 1:ord
                for j = 1:ord
                    fig_idx=j+fig;
                    data = squeeze(DjPj_(i*j, :));
                    approxData = squeeze(Happrox(i*j, :));
                    figure(fig_idx); clf;
                    for subp = 1:2
                        if subp == 1
                            subplot(2, 1, subp)
                            semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real DjPj_{%d,%d}, mode %d', i,j, o)); hold all;
                            semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpoles_thismode)); hold all;
                            legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                        else
                            subplot(2, 1, subp)
                            semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag DjPj_{%d,%d}, mode %d ', i,j, o)); hold all;
                            semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpoles_thismode)); hold all;
                            legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                        end
                    end
                end
                fig=fig_idx;
            end
        end
    end
end

%%  Write the poles and residues of Yc and H in .dat

% Poles of Yc
polYc = fit_data_Yc.A;
polYc_dupes=find(imag(polYc)<0);
polYc(polYc_dupes)=[];
polYc_dupes = sort(polYc_dupes,'descend');

% Residues of Yc
if polYc_dupes~=0
    for k=1:length(polYc_dupes)
        SERYc.C(:,polYc_dupes(k))=[];
    end
end

% Poles of H
for i=1:size(poles_wrapper,2)
    Poles_Pj(i,:) = poles_wrapper(:,i);
    Poles_Pj_dupes{i}=find(imag(Poles_Pj{i,1})<0);
    Poles_Pj{i,1}(Poles_Pj_dupes{1,i},:) = [];
end

% Residues of Yc
for k=1:length(Poles_Pj_dupes)
    polD=Poles_Pj_dupes{1,k}~=0;
    if polD
        rDjPj{1,k}(:,Poles_Pj_dupes{1,k})=[];
    end
end

npYc=size(polYc,1);

filename = fullfile(currPath,[jobid '_fitULM.dat']);

fid = fopen(filename,'wt');
fprintf(fid,'%d\n',ord); 	%Number of phases
fprintf(fid,'%d\n',ord); 	%Number of de modes
fprintf(fid,'%d\n',npYc); 	%Number of poles of Yc

for k=1:ord
    fprintf(fid,'%d\n',size(Poles_Pj{k,1},1));
end

%Optimal time delay
for jj = 1:size(Optimal_time,1)
    fprintf(fid,'%.16E\n',Optimal_time(jj));
end

% Write the poles of Yc and their residues
for k=1:npYc
    % Write the pole
    pol=polYc(k);
    if isreal(pol)
        if pol>0
            fprintf(fid,'%.16E\n',-(pol)); 	% Forces the real pole to be negative
        else
            fprintf(fid,'%.16E\n',pol);
        end
    else
        fprintf(fid,'+%.16E\n%.16E\n',abs(real(pol)),imag(pol));	% Forces the complex pole to have a positive real part
    end

    % Write the residues for the corresponding pole
    res=SERYc.C(:,k);
    res=reshape(res, [ord ord]);
    res=res.';
    res=triu(res);
    res=reshape(res.', [ord*ord 1]);
    res(res==0)=[];
    nres=length(res);
    for kk=1:nres
        if isreal(pol)
            fprintf(fid,'%.16E\n',real(res(kk)));
        else
            fprintf(fid,'%.16E\n%.16E\n',real(res(kk)),imag(res(kk)));
        end
    end
end

% Write the poles of DjPj and their residues for each mode
for k=1:ord
    for j=1:size(Poles_Pj{k,1},1)
        pol=Poles_Pj{k,1}(j);
        if isreal(pol)
            if pol>0
                fprintf(fid,'%.16E\n',-(pol));
            else
                fprintf(fid,'%.16E\n',pol);
            end
        else
            fprintf(fid,'+%.16E\n%.16E\n',abs(real(pol)),imag(pol));
        end
        res = rDjPj{1,k}(:,j);
        nres=length(res);
        for kk=1:nres
            if isreal(pol)
                fprintf(fid,'%.16E\n',real(res(kk)));
            else
                fprintf(fid,'%.16E\n%.16E\n',real(res(kk)),imag(res(kk)));
            end
        end

    end
end

% Independent terms of Yc
k0=SERYc.D;
k0=reshape(k0, [ord ord]);
k0=k0.';
k0=triu(k0);
k0=reshape(k0.', [ord*ord 1]);
k0(k0==0)=[];
fprintf(fid,'%.16E\n',real(k0));

fclose(fid);
end % of main function

function [H_mod]=calc_prop_fun(Ti,g,length,ord,num_files)
Ti_dis=zeros(ord,ord);
H_mod=zeros(num_files,ord);
for k=1:1:num_files
    for o=1:ord
        Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
    end

    H_mod(k,:)=exp(-g(k,:)*length);

end
end % of auxiliary propagation function

function [tau_opt, Aj_shifted] = findoptimtau(f,vel,H,line_length,optsRW)
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
    tmin = 0.95*tmin;
end

options = optimset('Display','none', 'TolX',1e-6);
tau_opt=fminbnd(@(x)fcalc(H,f,x,optsRW),tmin,tmax,options);
Aj_shifted = H.*exp(1i.*2.*pi.*f.*tau_opt);
end % of auxiliary function find optimum tau

function [out] = fcalc(H, f, tau,optsRW)
fun = H.*exp(1i.*2.*pi.*f.*tau);
[~, ~, ~, err_dB] = rationalfit_wrapper(fun,f,optsRW);
out =  err_dB;
end % yo dawg, I made a mess inside a mess so you can mess while u mess

%backup
% if ~isfield(optsRW,'plot1');optsRW.plot1=false;end;
% if ~isfield(optsRW,'plot2');optsRW.plot2=false;end;
% if ~isfield(optsRW,'Niter');optsRW.Niter = [20 50];end;
% if ~isfield(optsRW,'TendsToZero');optsRW.TendsToZero = true;end;
% if ~isfield(optsRW,'Npoles');optsRW.Npoles = [2 20];end;
% 
% % Fit_print= false;    %Flag for plot for all fits
% transp=@(x) x.';
% 
% %% Initial data
% 
% if size(f,1) == 1
%     frequency=transp(f);
% else
%     frequency=f;
% end
% 
% w = 2*pi.*frequency;    % Calculation of angular frequency (rad/sec)
% s = transp(1j.*w);      % Calculation of the complex angular frequency
% freq_siz = size(frequency,1);   % Number of frequency samples
% ord=size(Z,2);  % Order of the system
% 
% %% Calculating the frequency-dependent transformation matrix
% 
% % Raw impedance and admittance matrices
% if size(Z,3)==freq_siz
%     Z=permute(Z,[3 1 2]);
%     Y=permute(Y,[3 1 2]);
% end
% 
% [modif_T,g_dis]=LM_calc_norm_str(ord,freq_siz,Z,Y,frequency); % Transformation matrix using the Levenberg-Marquardt (LM)
% 
% Ti = zeros(size(Z,3), size(Z,3), size(frequency,1));
% invTi = zeros(size(Z,3), size(Z,3), size(frequency,1));
% for k = 1:freq_siz
%     for o = 1:ord
%         Ti(o,:,k) = modif_T(k,(o-1)*ord+1:o*ord); % Ti in N_phases x N_phases x N_samples
%     end
%     invTi(:,:,k) = inv(Ti(:,:,k)); % Calculation the inverse of Ti
% end
% 
% %% Calculating the characteristic impedance, admittance and the trace your trace
% 
% [~,~,~,Ych] = calc_char_imped_admit(modif_T,Z,Y,ord,freq_siz); % Calculation of characteristic impedance and admittance
% 
% Yc = zeros(size(Z,3), size(Z,3), size(frequency,1));
% for k = 1:freq_siz
%     for o = 1:ord
%         Yc(o,:,k) = Ych(k,(o-1)*ord+1:o*ord); %  Ych in N_phases x N_phases x N_phases
%     end
% end
% 
% traceYc = zeros(size(frequency));
% Yc_t = zeros(size(Z,3), size(Z,3), size(frequency,1));
% for k = 1 : size(Z,1)
%     Yc_t(:,:,k) = transp(Yc(:,:,k));
%     traceYc(k) = trace( Yc(:,:,k) ); % Calculating the Yc trace
% end
% 
% %% Calculating the propagation function and phase velocity
% [Aj] = calc_prop_fun(modif_T,g_dis,line_length,ord,freq_siz); 	% Propagation function for each modes
% 
% vel = zeros(size(frequency,1), size(Z,3));
% for m = 1:ord
%     vel(:,m) = (2*pi*frequency)./imag(g_dis(:,m)); 	% Calculation of phase velocity
% end
% %% Calculating Yc trace poles
% fitYc = transp(traceYc);
% [fit_data_Yc, ffit, numpolYc, ~] = rationalfit_wrapper(fitYc,frequency,optsRW); %Calculating the poles
% 
% if ~ispassive(fit_data_Yc)
%     warning('Resulting fitted traceYc is not passive. Results might be numerically unstable or inaccurate.')
% end
% 
% if optsRW.plot1
%     fig=1;
%     figure(fig);
%     semilogx(frequency,abs(traceYc),'-b', 'linewidth', 3); hold on
%     semilogx(frequency,abs(ffit),'--r', 'linewidth', 3)
%     legend('trace Yc', 'Fit trace Yc', 'location', 'north' )
%     ylabel('Magnitude')
%     xlabel('Frequency [Hz]')
%     grid on
% end
% %% Calculating the residuals of Yc
% 
% % Opts for Vector Fitting
% clear opts
% opts.relax = 1;
% opts.stable = 1;
% opts.asymp = 2;
% opts.skip_pole = 1;
% opts.skip_res = 0;
% opts.cmplx_ss = 1;
% opts.spy1 = 0;
% opts.spy2 = 0;
% 
% rYc = reshape(Yc_t, size(Yc_t,2)*size(Yc_t,1), size(Yc_t,3)); % Yc_t in  N_phases.N_phases x N_samples
% 
% for j = 1:1 % número de iterações
%     weight = ones(1,numel(s)); % peso - 1 x Namostras
%     [SERYc,~,~,Ycapprox,~] = vectfit3(rYc, s, transp(fit_data_Yc.A), weight, opts); % Function VF by Gustavsen
% end
% 
% Ycapprox_ = zeros(size(Z,3), size(Z,3), size(frequency,1));
% for k = 1:freq_siz
%     for o = 1:ord
%         Ycapprox_(o,:,k) = Ycapprox((o-1)*ord+1:o*ord,k); % Ycapprox in N_phases x N_phases x N_samples
%     end
% end
% 
% % Plot for each element of the Yc matrix (N_phases x N_phases)
% if optsRW.plot1
%     linewidth = 3;
%     approxColor = ':r';
%     for i = 1:ord
%         for j = 1:ord
%             fig_idx=fig+j;
%             data = squeeze(Ych(:,i*j));
%             approxData = squeeze(Ycapprox(i*j, :));
%             figure(fig_idx); clf;
%             for subp = 1:2
%                 if subp == 1
%                     subplot(2, 1, subp)
%                     semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Yc_{%d,%d}', i, j)); hold all;
%                     semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpolYc)); hold all;
%                     legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
% 
%                 else
%                     subplot(2, 1, subp)
%                     semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Yc_{%d,%d}', i, j)); hold all;
%                     semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpolYc)); hold all;
%                     legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
%                 end
%             end
%         end
%         fig=fig_idx;
%     end
% end
% %% Calculating Pj with optimal time delay
% 
% for m = 1:ord
%     [tau_opt, Aj_shifted] = findoptimtau(frequency,vel(:,m),Aj(:,m),line_length,optsRW); 	% Optimized time finding function
%     Optimal_time(m,1) = tau_opt;
%     Pj(m,:)= Aj_shifted;
% end
% 
% P_j_t = zeros(length(frequency),ord);
% for i = 1:ord
%     for j = 1:length(frequency)
%         P_j_t(j,i) = transp(Pj(i,j));
%     end
% end
% 
% Pjj = zeros(size(P_j_t,2),size(P_j_t,2),size(f,2));
% for k=1:freq_siz
%     Pjj(:,:,k)=diag(P_j_t(k,:)); % Pj in a diagonal matrix
% end
% 
% %% Calculating fitting of Pj
% 
% for t=1:size(Pj,1)
%     fitPj=transp(Pj(t,:));
%     [fit_data, ffit, ~, ~] = rationalfit_wrapper(fitPj,frequency,optsRW);	%Calculating the poles for each mode
%     poles_wrapper(t) = {fit_data.A};
%     Pjapprox(t,:) = ffit;
% end
% 
% %	Plot for the Pj adjustment poles
% if optsRW.plot1
%     fig=1;
%     linewidth = 3;
%     approxColor = ':r';
%     for i = 1:ord
%         fig_idx=i+fig;
%         data = Pj(i,:);
%         approxData = Pjapprox(i,:);
%         figure(fig_idx); clf;
%         for subp = 1:2
%             if subp == 1
%                 subplot(2, 1, subp)
%                 semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Pj_{mode %d}', i)); hold all;
%                 semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', size(poles_wrapper{1},1))); hold all;
%                 legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
% 
%             else
%                 subplot(2, 1, subp)
%                 semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Pj_{mode %d}', i)); hold all;
%                 semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', size(poles_wrapper{1},1))); hold all;
%                 legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
%             end
%         end
%     end
%     fig=fig_idx;
% end
% 
% %% Calculating the residues of DjPj
% for m=1:ord
%     for k=1:freq_siz
%         mode(m).D(:,:,k)=Ti(:,m,k)*invTi(m,:,k); % Idempotent Dj matrix -> Eliminates indeterminacy problems
%     end
% end
% 
% %	Opts for Vector fitting
% clear opts
% opts.relax = 1;
% opts.stable = 1;
% opts.asymp = 1;
% opts.skip_pole = 1;
% opts.skip_res = 0;
% opts.cmplx_ss = 1;
% opts.spy1 = 0;
% opts.spy2 = 0;
% opts.passive = 0;
% 
% for o = 1:ord
%     clear DjPj
%     for k = 1:freq_siz
%         DjPj(:,:,k) = mode(o).D(:,:,k)*Pjj(o,o,k);
%         DjPj_t(:,:,k) = transp(DjPj(:,:,k));
%     end
%     DjPj_ = reshape(DjPj_t, size(DjPj_t,2)*size(DjPj_t,1), size(DjPj_t,3)); %Yc_t in (N_phases*N_phases) x N_samples
%     poles = poles_wrapper{o};
%     numpoles_thismode = length(poles);
%     for j = 1:1
%         weight = ones(1,numel(s));
%         [SER_DjPj,~,~,Happrox,~] = vectfit3(DjPj_, s, poles, weight, opts); %
%         rDjPj(o) = {SER_DjPj.C};	% Residues for DjPj
% 
%         if optsRW.plot1
%             linewidth = 3;
%             approxColor = ':r';
%             for i = 1:ord
%                 for j = 1:ord
%                     fig_idx=j+fig;
%                     data = squeeze(DjPj_(i*j, :));
%                     approxData = squeeze(Happrox(i*j, :));
%                     figure(fig_idx); clf;
%                     for subp = 1:2
%                         if subp == 1
%                             subplot(2, 1, subp)
%                             semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real DjPj_{%d,%d}, mode %d', i,j, o)); hold all;
%                             semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpoles_thismode)); hold all;
%                             legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
% 
%                         else
%                             subplot(2, 1, subp)
%                             semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag DjPj_{%d,%d}, mode %d ', i,j, o)); hold all;
%                             semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', numpoles_thismode)); hold all;
%                             legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;
% 
%                         end
%                     end
%                 end
%                 fig=fig_idx;
%             end
%         end
%     end
% end
% 
% %%  Write the poles and residues of Yc and H in .dat
% 
% % Poles of Yc
% polYc = fit_data_Yc.A;
% polYc_dupes=find(imag(polYc)<0);
% polYc(polYc_dupes)=[];
% 
% % Residues of Yc
% if polYc_dupes~=0
%     for k=1:length(polYc_dupes)
%         SERYc.C(:,polYc_dupes(k))=[];
%     end
% end
% 
% % Poles of H
% for i=1:size(poles_wrapper,2)
%     Poles_Pj(i,:) = poles_wrapper(:,i);
%     Poles_Pj_dupes{i}=find(imag(Poles_Pj{i,1})<0);
%     Poles_Pj{i,1}(Poles_Pj_dupes{1,i},:) = [];
% end
% 
% % Residues of Yc
% for k=1:length(Poles_Pj_dupes)
%     polD=Poles_Pj_dupes{1,k}~=0;
%     if polD
%         rDjPj{1,k}(:,Poles_Pj_dupes{1,k})=[];
%     end
% end
% 
% npYc=size(polYc,1);
% 
% filename = fullfile(currPath,[jobid '_fitULM.dat']);
% 
% fid = fopen(filename,'wt');
% fprintf(fid,'%d\n',ord); 	%Number of phases
% fprintf(fid,'%d\n',ord); 	%Number of de modes
% fprintf(fid,'%d\n',npYc); 	%Number of poles of Yc
% 
% for k=1:ord
%     fprintf(fid,'%d\n',size(Poles_Pj{k,1},1));
% end
% 
% %Optimal time delay
% for jj = 1:size(Optimal_time,1)
%     fprintf(fid,'%.16E\n',Optimal_time(jj));
% end
% 
% % Write the poles of Yc and their residues
% for k=1:npYc
%     % Write the pole
%     pol=polYc(k);
%     if isreal(pol)
%         if pol>0
%             fprintf(fid,'%.16E\n',-(pol)); 	% Forces the real pole to be negative
%         else
%             fprintf(fid,'%.16E\n',pol);
%         end
%     else
%         fprintf(fid,'+%.16E\n%.16E\n',abs(real(pol)),imag(pol));	% Forces the complex pole to have a positive real part
%     end
% 
%     % Write the residues for the corresponding pole
%     res=SERYc.C(:,k);
%     res=reshape(res, [ord ord]);
%     res=res.';
%     res=triu(res);
%     res=reshape(res.', [ord*ord 1]);
%     res(res==0)=[];
%     nres=length(res);
%     for kk=1:nres
%         if isreal(pol)
%             fprintf(fid,'%.16E\n',real(res(kk)));
%         else
%             fprintf(fid,'%.16E\n%.16E\n',real(res(kk)),imag(res(kk)));
%         end
%     end
% end
% 
% % Write the poles of DjPj and their residues for each mode
% for k=1:ord
%     for j=1:size(Poles_Pj{k,1},1)
%         pol=Poles_Pj{k,1}(j);
%         if isreal(pol)
%             if pol>0
%                 fprintf(fid,'%.16E\n',-(pol));
%             else
%                 fprintf(fid,'%.16E\n',pol);
%             end
%         else
%             fprintf(fid,'+%.16E\n%.16E\n',abs(real(pol)),imag(pol));
%         end
%         res = rDjPj{1,k}(:,j);
%         nres=length(res);
%         for kk=1:nres
%             if isreal(pol)
%                 fprintf(fid,'%.16E\n',real(res(kk)));
%             else
%                 fprintf(fid,'%.16E\n%.16E\n',real(res(kk)),imag(res(kk)));
%             end
%         end
% 
%     end
% end
% 
% % Independent terms of Yc
% k0=SERYc.D;
% k0=reshape(k0, [ord ord]);
% k0=k0.';
% k0=triu(k0);
% k0=reshape(k0.', [ord*ord 1]);
% k0(k0==0)=[];
% fprintf(fid,'%.16E\n',real(k0));
% 
% fclose(fid);
% end
% 
% function [H_mod]=calc_prop_fun(Ti,g,length,ord,num_files)
% 
% Ti_dis=zeros(ord,ord);
% H_mod=zeros(num_files,ord);
% 
% for k=1:1:num_files
% 
%     for o=1:ord
%         Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
%     end
% 
%     H_mod(k,:)=exp(-g(k,:)*length);
% 
% end
% end
% 
% function [tau_opt, Aj_shifted] = findoptimtau(f,vel,H,line_length,optsRW)
% j=length(f)-1;
% w=2*pi*f;
% Ns=length(w);
% absH=abs(H);
% 
% phase1=(pi/2)*log((absH(j+1)/absH(j-1)))/(log(w(j+1)/w(j-1)));
% phase2=0;
% term2=log((absH(j+1)/absH(j-1))) /(log(w(j+1)/w(j-1)));
% for k=2:Ns-1
%     term1=log(absH(k+1)/absH(k-1)) /(log(w(k+1)/w(k-1)));
%     if k~=j
%         phase2=phase2+(abs(term1)-abs(term2))*log(coth(abs(log(w(k)/w(j)))/2))*log(w(k+1)/w(k));
%     end
% end
% 
% phase2=phase2/pi;
% ph=(phase1-phase2);
% tau_mps=(line_length/vel(j)) + ph/(w(j));
% 
% tmin=1*abs(tau_mps);
% tmax=1*abs(line_length/vel(end));
% 
% if tmin >= tmax
%     tmin = 0.95*tmin;
% end
% 
% options = optimset('Display','none', 'TolX',1e-6);
% tau_opt=fminbnd(@(x)fcalc(H,f,x,optsRW),tmin,tmax,options);
% Aj_shifted = H.*exp(1i.*2.*pi.*f.*tau_opt);
% 
% end
% 
% 
% function [out] = fcalc(H, f, tau,optsRW)
% fun = H.*exp(1i.*2.*pi.*f.*tau);
% [~, ~, ~, err_dB] = rationalfit_wrapper(fun,f,optsRW);
% out =  err_dB;
% end
