function [] = makeULModel(Z,Y,f,line_length,ERR,Npoles,jobid,currPath)

ZYprnt=true;

transp=@(x) x.';

%% Dados

if size(f,1) == 1
    frequency=f.';
else
    frequency=f;
end

w = 2*pi.*frequency; % % cálculo da frequência angular (rad/seg)
s = 1j.*w; % cálculo da frequência angular complexa
freq_siz = size(frequency,1); % número de amostras da frequência
ord=size(Z,2);

%% Calculando matriz de transformação dependente da frequência

% Raw impedance and admittance matrices
if size(Z,3)==freq_siz
    Z=permute(Z,[3 1 2]);
    Y=permute(Y,[3 1 2]);
end

[modif_T,g_dis]=LM_calc_norm_str(ord,freq_siz,Z,Y,frequency); % obtem as matrizes de transformação pelo Levenberg-Marquardt (LM)

Ti = zeros(size(Z,3), size(Z,3), size(frequency,1));
invTi = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1:freq_siz
    for o = 1:ord
        Ti(o,:,k) = modif_T(k,(o-1)*ord+1:o*ord); % reorganiza as dimensões de Ti para que seja Nfases x Nfases x Namostras
    end
    invTi(:,:,k) = inv(Ti(:,:,k)); % calcula a inversa de Ti
end

%% Calculando a impedância e admitância característica e o traço de Yc

[~,~,Zch,Ych] = calc_char_imped_admit(modif_T,Z,Y,ord,freq_siz); %calculo da impedancia e admitancia caracteristica

Yc = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1:freq_siz
    for o = 1:ord
        Yc(o,:,k) = Ych(k,(o-1)*ord+1:o*ord); % reorganiza as dimensões de Ych para que seja Nfases x Nfases x Namostras
    end
end

traceYc = zeros(size(frequency));
Yc_t = zeros(size(Z,3), size(Z,3), size(frequency,1));
for k = 1 : size(Z,1)
    Yc_t(:,:,k) = transp(Yc(:,:,k)); % transposta da dmitância característica para usar depois para que os resultados sejam no formato [11; 12; 13; 14; 15; 16; 21; 22; 23; 24; 25; 26; 31; 32; 33; 34; 35; 36 ... 61; 62; 63; 64; 65; 66] vetor coluna
    traceYc(k) = trace( Yc(:,:,k) ); % cálculo do traço de Yc - na dissertação do Zanon fala que os polos da admitância característica são calculados pelo seu traço (soma dos elementos da diagonal principal)
end

%% Calculando a função de propagação e a velocidade de fase
[Aj] = calc_prop_fun(modif_T,g_dis,line_length,ord,freq_siz); % função que calcula a função de propagação

vel = zeros(size(frequency,1), size(Z,3));
for m = 1:ord
    vel(:,m) = (2*pi*frequency)./imag(g_dis(:,m)); % calcula a velocidade de fase
end

%% Calculando polos do traço de Yc

% Npoles = numpoles; % número de polos
Ycapprox = cell(1, numel(Npoles));

% opts usadas no metodo do VF
opts.relax = 1;      % use ajuste vetorial com restrição de não trivialidade relaxada
opts.stable = 1;     % aplicar pólos estáveis
opts.asymp = 2;      % fitting com D~=0, E=0
opts.skip_pole = 0;  % não pule o cálculo dos polos
opts.skip_res = 0;   % não pule o cálculo dos resíduos
opts.cmplx_ss = 1;   % crie um modelo de espaço de estado complexo
opts.spy1 = 0;       % sem plotagem para o primeiro estágio do VF
opts.spy2 = 0;       % criar gráfico de magnitude para ajuste de f(s)

f = traceYc.'; % renomeando o traço de Yc e realizando sua transposta - 1 x Namostras
s = s(:).'; % realizando a transposta de s - 1 x Namostras

% for k = 1 : numel(Npoles)
% 
%     poles = linspace( frequency(1), frequency(end), Npoles(k) ); % polos iniciais - 1 x Namostras
% 
%     for j = 1:10 % número de iterações
%         weight = ones(1,numel(s)); % peso - 1 x Namostras
%         [SERtrYc,poles,~,fit,~] = vectfit3(f, s, poles, weight, opts); % função do VF elaborada por Gustavsen
%     end
%     fittedPoles = poles; % polos resultantes de Yc
% end

fit_data = rationalfit(frequency.',f,'NPoles',Npoles,'TendsToZero',false);
[ffit,freq] = freqresp(fit_data,frequency);
fittedPoles=fit_data.A.';

if ~ispassive(fit_data)
    warning('Resulting fitted traceYc is not passive. Results might be numerically unstable or inaccurate.')
end

if ZYprnt
    fig=1;
    figure(fig);
    semilogx(frequency,abs(traceYc),'-b', 'linewidth', 3); hold on
    semilogx(frequency,abs(ffit),'--r', 'linewidth', 3)
    legend('trace Yc', 'Fit trace Yc', 'location', 'north' )
    ylabel('Magnitude')
    xlabel('Frequency [Hz]')
    grid on
end
%% Calculando os resíduos de Yc

clear opts
opts.relax = 1;      % use ajuste vetorial com restrição de não trivialidade relaxada
opts.stable = 1;     % aplicar pólos estáveis
opts.asymp = 2;      % fitting com D~=0, E=0
opts.skip_pole = 1;  % pule o cálculo dos polos
opts.skip_res = 0;   % não pule o cálculo dos resíduos
opts.cmplx_ss = 1;   % crie um modelo de espaço de estado complexo
opts.spy1 = 0;       % sem plotagem para o primeiro estágio do VF
opts.spy2 = 0;       % criar gráfico de magnitude para ajuste de f(s)

for k = 1 : numel(Npoles)
    rYc = reshape(Yc_t, size(Yc_t,2)*size(Yc_t,1), size(Yc_t,3)); % reorganiza as dimensões de Yc_t para que fique Nfases.Nfases x Namostras % [11; 12; 13; 14; 15; 16; 21; 22; 23; 24; 25; 26; 31; 32; 33; 34; 35; 36 ... 61; 62; 63; 64; 65; 66] vetor coluna
    s = s(:).'; % realizando a transposta de s - 1 x Namostras

    poles = fittedPoles; % polos iniciais obtidos no traço de Yc

    for j = 1:1 % número de iterações
        weight = ones(1,numel(s)); % peso - 1 x Namostras
        [SERYc,poles,rmserr,Ycapprox,~] = vectfit3(rYc, s, poles, weight, opts); % função do VF elaborada por Gustavsen
    end
end

if ZYprnt

    linewidth = 3;
    approxColor = ':r';
    for i = 1:ord
        for j = 1:ord
            fig_idx=i*j+fig;
            data = squeeze(rYc(i*j, :));
            approxData = squeeze(Ycapprox(i*j, :));
            figure(fig_idx); clf;
            for subp = 1:2
                if subp == 1
                    subplot(2, 1, subp)
                    semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Yc_{%d,%d}', i, j)); hold all;
                    semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                else
                    subplot(2, 1, subp)
                    semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Yc_{%d,%d}', i, j)); hold all;
                    semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                    legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                end
            end
        end
    end
fig=fig+ord*ord;

end

%% Calculando Pj com tau otimizado

% opts usadas no metodo do VF
clear opts
opts.relax = 1;      % use ajuste vetorial com restrição de não trivialidade relaxada
opts.stable = 1;     % aplicar pólos estáveis
opts.skip_pole = 0;  % não pule o cálculo dos polos
opts.skip_res = 0;   % não pule o cálculo dos resíduos
opts.cmplx_ss = 1;   % crie um modelo de espaço de estado complexo
opts.spy1 = 0;       % sem plotagem para o primeiro estágio do VF
opts.spy2 = 0;       % criar gráfico de magnitude para ajuste de f(s)
opts.logx = 1;       % usar eixo logarítmico de abcissas
opts.logy = 1;       % usar eixo de ordenadas logarítmicas
opts.errplot = 1;    % incluir desvio no gráfico de magnitude
opts.phaseplot = 0;  % também produz plotagem do ângulo de fase (além da magnitude)
opts.legend = 1;     % inclua legendas nos gráficos
opts.asymp = 1;      % fitting com D=0, E=0
opts.weightscheme = 3; % pesos 1 = uniform weights; 2 = inv(sqrt(norm)) 3 = inv(norm)
opts.firstguesstype = 1; % polos de valor real - usado no VF wrapper
opts.output_messages = false; % - usado no VF wrapper
opts.passive = 0; % - usado no VF wrapper

% ERR = .1/100; % tolerância do erro
for m = 1:ord
    [tau_opt, fun] = findoptimtau(frequency,vel(:,m),Aj(:,m),line_length,ERR,opts); % função que faz o VF da função de propagação considerando o tau ótimo, nesse caso uso para obter o Pj com tau otimizado
%     [~, ~, ~, ~, ~, ~, ~, tau_opt, ~, fun] = findoptimtau(frequency,vel(:,m),Aj(:,m),line_length,ERR,opts); % função que faz o VF da função de propagação considerando o tau ótimo, nesse caso uso para obter o Pj com tau otimizado

    fitOHLT_H(m).tau_opt = tau_opt; % tau otimizado
    fitOHLT_H(m).fun = fun; % Pj = H.*exp(1i.*2.*pi.*f.*tau) -  eq. 3.22 dissertação zanon;
end

P_j = zeros(ord,length(frequency));
P_j_t = zeros(length(frequency),ord);
for i = 1:ord
    for j = 1:length(frequency)
        P_j(i,j) = fitOHLT_H(i).fun(j,:); % passando Pj para forma Nfases x Namostras para usar no fitting
        P_j_t(j,i) = P_j(i,j).'; % transposta de Pj
    end
end

Pj = zeros(size(Yc));
for k=1:freq_siz
    Pj(:,:,k)=diag(P_j_t(k,:)); % transformando Pj em uma matriz diagonal para usar depois no cálculo de DjPj
end

%% Calculando fitting de Pj

clear opts
opts.relax = 1;      % use ajuste vetorial com restrição de não trivialidade relaxada
opts.stable = 1;     % aplicar pólos estáveis
opts.asymp = 1;      % fitting com D=0, E=0
opts.skip_pole = 0;  % nao pule o cálculo dos polos
opts.skip_res = 0;   % nao pule o cálculo dos resíduos
opts.cmplx_ss = 1;   % crie um modelo de espaço de estado complexo
opts.spy1 = 0;       % sem plotagem para o primeiro estágio do VF
opts.spy2 = 0;       % criar gráfico de magnitude para ajuste de f(s)

% for k = 1 : numel(Npoles)
%     s = s(:).'; % realizando a transposta de s - 1 x Namostras
%     Pjj = P_j; % renomeando P_j Nfases x Namostras
% 
%     poles = linspace( frequency(1), frequency(end), Npoles(k) ); % polos iniciais - 1 x Namostras
% 
%     for j = 1:10 % número de iterações
%         weight = ones(1,numel(s)); % peso - 1 x Namostras
%         [SERPj,poles,~,Pjapprox,~] = vectfit3(Pjj, s, poles, weight, opts); % função do VF elaborada por Gustavsen
%     end
%     polesPj = poles; % polos resultantes de Yc
% end

fit_data = rationalfit(frequency.',P_j.','NPoles',Npoles,'TendsToZero',true);
foundpassive=false;
for i=1:length(fit_data)
    [ffit,freq] = freqresp(fit_data(i),frequency);
    Pjapprox(i,:)=ffit;
    polesPj=fit_data(i).A.';

    if ispassive(fit_data(i))
        passivemode=i;
        foundpassive=true;
    else
        warning('Resulting P_j fit for mode %d is not passive. Results might become a dumpster fire.',i)
        foundpassive=false;
    end
    
    if foundpassive;polesPj=fit_data(passivemode).A.';end

end




if ZYprnt

    linewidth = 3;
    approxColor = ':r';
    for i = 1:ord
        fig_idx=i+fig;
        data = P_j(i,:);
        approxData = Pjapprox(i,:);
        figure(fig_idx); clf;
        for subp = 1:2
            if subp == 1
                subplot(2, 1, subp)
                semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real Pj_{mode %d}', i)); hold all;
                semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

            else
                subplot(2, 1, subp)
                semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag Pj_{mode %d}', i)); hold all;
                semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

            end
        end
    end
fig=fig+ord;

end

%% Calculando os resíduos de DjPj

for m=1:ord
    for k=1:freq_siz
        mode(m).D(:,:,k)=Ti(:,m,k)*invTi(m,:,k); % matriz Dj idempotentes
    end
end

clear opts
opts.relax = 1;      % use ajuste vetorial com restrição de não trivialidade relaxada
opts.stable = 1;     % aplicar pólos estáveis
opts.asymp = 1;      % fitting com D=0, E=0
opts.skip_pole = 1;  % pule o cálculo dos polos
opts.skip_res = 0;   % não pule o cálculo dos resíduos
opts.cmplx_ss = 1;   % crie um modelo de espaço de estado complexo
opts.spy1 = 0;       % sem plotagem para o primeiro estágio do VF
opts.spy2 = 0;       % criar gráfico de magnitude para ajuste de f(s)

for o = 1:ord
    clear DjPj
    for k = 1:freq_siz
        DjPj(:,:,k) = mode(o).D(:,:,k)*Pj(o,o,k); % Dj.Pj - eq. 3.23 dissertação zanon
        DjPj_t(:,:,k) = transp(DjPj(:,:,k)); % transposta de DjPj
    end
    for k = 1 : numel(Npoles)
        rDjPj = reshape(DjPj_t, size(DjPj_t,2)*size(DjPj_t,1), size(DjPj_t,3)); % reorganiza as dimensões de Yc_t para que fique (Nfases*Nfases) x Namostras % [11; 12; 13; 14; 15; 16; 21; 22; 23; 24; 25; 26; 31; 32; 33; 34; 35; 36 ... 61; 62; 63; 64; 65; 66] vetor coluna
        s = s(:).'; % realizando a transposta de s - 1 x Namostras

        poles = polesPj; % polos iniciais obtidos no Pj

        for j = 1:1
            weight = ones(1,numel(s)); % peso - 1 x Namostras
            [SERH,poles,~,Happrox,~] = vectfit3(rDjPj, s, poles, weight, opts); % função do VF elaborada por Gustavsen

            Cij(:,:,o) = SERH.C;

            if ZYprnt

                linewidth = 3;
                approxColor = ':r';
                for i = 1:ord
                    for j = 1:ord
                        fig_idx=o*i*j+fig;
                        data = squeeze(rDjPj(i*j, :));
                        approxData = squeeze(Happrox(i*j, :));
                        figure(fig_idx); clf;
                        for subp = 1:2
                            if subp == 1
                                subplot(2, 1, subp)
                                semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real DjPj_{%d}, mode %d', i*j, o)); hold all;
                                semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                                legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                            else
                                subplot(2, 1, subp)
                                semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag DjPj_{%d}, mode %d', i*j, o)); hold all;
                                semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                                legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                            end
                        end
                    end
                end

            end

        end
    end
end

%% Para escrever os polos e residuos de Yc e H em txt

ii = 1;
for m = 1:ord
    [b, ~] = size(fitOHLT_H(m).tau_opt); % reorganiza as dimensões de fitOHLT_H(m).tau_opt para que fique número de taus total (somando todos os modos) x 1 - ex LT trifásica (3 x 1)
    D_aux(:,1) =  fitOHLT_H(m).tau_opt;
    for i = 1:b
        D(ii,1) = D_aux(i,1);
        ii = ii + 1;
    end
    clear D_aux
end

pol = fittedPoles.';
polPj = polesPj.';
rYc = reshape(SERYc.C, size(SERYc.C,1)*size(SERYc.C,2), 1); % reorganiza as dimensões de SERYc.C para que fique nfases.nfases.npolos x 1
rCij = reshape(Cij, size(Cij,1)*size(Cij,2), size(Cij,3)); % reorganiza as dimensões de Cij para que fique nfases.nfases.npolos x nmodos
rresiduos = reshape(rCij, size(rCij,1)*size(rCij,2), 1); % reorganiza as dimensões de SERYc.C para que fique nfases.nfases.npolos.modos x 1

filename = fullfile(currPath,[jobid '_fitULM2.dat']);

fid = fopen(filename,'wt');
fprintf(fid,'%d\n',size(Yc,1)); %numero de fases
fprintf(fid,'%d\n',size(Yc,1)); %numero de modos
fprintf(fid,'%d\n',size(pol,1)); %numero de polos Yc
fprintf(fid,'%d\n',size(polPj,1)); %numero de polos A

%polYcvet
for ii = 1:size(pol,1)
    fprintf(fid,'%.16e\t%.16e\n',real(pol(ii)),imag(pol(ii)));
end
%resYc
for kk = 1:size(rYc,1)
    fprintf(fid,'%.16e\t%.16e\n',real(rYc(kk)),imag(rYc(kk)));
end
%polA
for ii = 1:size(polPj,1)
    fprintf(fid,'%.16e\t%.16e\n',real(polPj(ii)),imag(polPj(ii)));
end
%resA
for kk = 1:size(rresiduos,1)
    fprintf(fid,'%.16e\t%.16e\n',real(rresiduos(kk)),imag(rresiduos(kk)));
end
% taus otimizados
for jj = 1:size(D,1)
    fprintf(fid,'%.16e\n',D(jj));
end
% resid0Yc
for ii = 1:(size(SERYc.D,1))
    fprintf(fid,'%.16e\t%.16e\n',real(SERYc.D(ii)),imag(SERYc.D(ii)));
end
fclose(fid);
end

function [H_mod]=calc_prop_fun(Ti,g,length,ord,num_files)

Ti_dis=zeros(ord,ord);
H_mod=zeros(num_files,ord);

for k=1:1:num_files

    for o=1:ord
        Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
    end

    H_mod(k,:)=exp(-g(k,:)*length);

end
end

% function [rmserr, pol, res, ks, NORD, ffit, tau_mps, tau_opt, tau, fun] = findoptimtau(f,vel,H,line_length,ERR,opts)
function [tau_opt, fun] = findoptimtau(f,vel,H,line_length,ERR,opts)
j=length(f)-1;

ff=f;

w=2*pi*ff;
Ns=length(w);

absH=abs(H);

w=2*pi*ff;

%First term in (9): Accurate_transmission_line_modeling_through_optima.pdf
phase1=(pi/2)*log((absH(j+1)/absH(j-1)))/(log(w(j+1)/w(j-1)));
%Second term in (9):
phase2=0;
term2=log((absH(j+1)/absH(j-1))) /(log(w(j+1)/w(j-1)));
for k=2:Ns-1
    term1=log(absH(k+1)/absH(k-1)) /(log(w(k+1)/w(k-1)));
    if k~=j
        phase2=phase2+(abs(term1)-abs(term2))*log(coth(abs(log(w(k)/w(j)))/2))*log(w(k+1)/w(k));
    end
end

phase2=phase2/pi;
ph=(phase1+phase2);
if ph<0
    ph=ph+(2*pi);
end
phase_min=ph; %Phase angle [rad]
tau_mps=(line_length/vel(j)) + phase_min/(w(j));


tau_a=line_length/3e8;
tau_b=(line_length/vel(j));


tmin=.8*abs(tau_a);
tmax=1.2*abs(tau_b);

% DEBUGME
% dt=abs(tmin-tmax)/100;
% ttau=tmin:dt:tmax
% % opts.NORD=20;
% % opts.Niter=20;
% for i=1:length(ttau)
%     [rmserr, pol, res, ks, NORD, ffit]=fcalc(H, f, ttau(i), ERR, opts);pol
%     err_tau(i)=rmserr;
% end
% figure;plot(ttau,err_tau)

options = optimset('Display','none', 'TolX',1e-6);
tau_opt=fminbnd(@(x)fcalc(H,ff,x,ERR,opts),tmin,tmax,options);
fun = H.*exp(1i.*2.*pi.*ff.*tau_opt);
% opts.NORD=20;
% opts.Niter=20;


% while 1

% tau_opt=fminbnd(@(x)fcalc(H,ff,x,ERR,opts),tmin,tmax,options);
%     tau_opt=tau_b;
%[rmserr, pol, res, ks, NORD, ffit, fun]=fcalc(H, ff, tau_opt, ERR, opts);

%     if max(imag(pol))  ~= 0
%         opts.NORD = opts.NORD-1;
%     else
%         break
%     end
% end

end


function [out] = fcalc(H, f, tau, ERR, opts)
fun = H.*exp(1i.*2.*pi.*f.*tau);
[fit_data,err] = rationalfit(f.',fun.','NPoles',20,'TendsToZero',true);

% [pol, res, ks, NORD, ffit, err] = vectfit_wrapper(fun,f,ERR,opts);
out=err;

end

% function [out, pol, res, ks, NORD, ffit, fun] = fcalc(H, f, tau, ERR, opts)
% fun = H.*exp(1i.*2.*pi.*f.*tau);
% [fit_data,err] = rationalfit(f.',fun.','NPoles',20,'TendsToZero',true);
% 
% % [pol, res, ks, NORD, ffit, err] = vectfit_wrapper(fun,f,ERR,opts);
% out=err;
% 
% end

