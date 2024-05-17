function [] = makeULModel3(Z,Y,f,line_length,ERR,Npoles,jobid,currPath)

ZYprnt= false;
group_modes = true;

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

if group_modes
    dtau_min=2*pi*10/360/frequency(end);
    % tau_modes_m=min(line_length./vel);
    tau_modes_m=line_length./vel(end,:);
    modes_list=[1:ord];
    num_modes=length(modes_list);
    k=1;
    while length(modes_list)>0
        this_tau=tau_modes_m(modes_list(1));
        dtau=abs(tau_modes_m-this_tau);
        same_group=find(dtau<dtau_min);
        Aj_r(:,k)=mean(Aj(:,same_group),2);
        vel_r(k,:)=mean(vel(:,same_group),2);
        % tau_r(:,k)=min(tau_modes_m(:,same_group)); %probably will not be used - debug purposes only
        modes_list(ismember(modes_list,same_group))=[];
        k=k+1;
    end
    ord=size(Aj_r,2);
    Aj=Aj_r;
    vel=vel_r.';
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

 fit_data_Yc = rationalfit(frequency.',f,'NPoles',Npoles,'TendsToZero',false);

[ffit,freq] = freqresp(fit_data_Yc,frequency);
fittedPoles=fit_data_Yc.A.';

if ~ispassive(fit_data_Yc)
    warning('Resulting fitted traceYc is not passive. Results might be numerically unstable or inaccurate.')
end

if ZYprnt
    fig=1;
    figure(fig);
    semilogx(frequency,abs(traceYc),'-b', 'linewidth', 3); hold on
    semilogx(freq,abs(ffit),'--r', 'linewidth', 3)
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
            fig_idx=fig+j;
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
        fig=fig_idx;
    end
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
opts.output_messages = false; % - usado no VF wrapper
opts.passive = 0; 

for m = 1:ord
    [tau_opt, fun] = findoptimtau(frequency,vel(:,m),Aj(:,m),line_length,ERR,opts,Npoles); % função que faz o VF da função de propagação considerando o tau ótimo, nesse caso uso para obter o Pj com tau otimizado
    D(m,1) = tau_opt;
    fun_(m,:)= fun;
end

% if group_modes
%     dtau_min=2*pi*10/360/frequency(end);
%     tau_modes_m=D;
%     modes_list=[1:ord].';
%     k=1;
%     while length(modes_list)>0
%         this_tau=tau_modes_m(modes_list(1));
%         dtau=abs(tau_modes_m-this_tau);
%         same_group=find(dtau<dtau_min);
%         fun_r(k,:)=mean(fun_(same_group,:),1);
%         vel_r(k,:)=mean(vel(:,same_group),2);
%         tau_r(k,:)=min(tau_modes_m(same_group,:));
%         modes_list(ismember(modes_list,same_group))=[];
%         k=k+1;
%     end
%     ord=size(fun_r,1);
%     fun_=fun_r;
%     D=tau_r;
% end


P_j = zeros(ord,length(frequency));
P_j_t = zeros(length(frequency),ord);
for i = 1:ord
    for j = 1:length(frequency)
        P_j(i,j) = fun_(i,j); 
        P_j_t(j,i) = P_j(i,j).'; % transposta de Pj
    end
end

Pj = zeros(size(P_j_t,2),size(P_j_t,2),size(f,2));
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


[fit_data_Pj,err_Pj] = rationalfit(frequency.',P_j.','NPoles',Npoles,'TendsToZero',true,'IterationLimit',[24 40],'DelayFactor',0.5);
if err_Pj>-40
       warning('Resulting Pj fit may be poor. Thread lightly.') 
end

for i=1:length(fit_data_Pj)
    [ffit,freq] = freqresp(fit_data_Pj(i),frequency);
    Pjapprox(i,:)=ffit;
    polesPj=fit_data_Pj(i).A.';
end

%Verifica a qualidade dos polos do ajuste de Pj
if ZYprnt
    fig=1;
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
 fig=fig_idx;
end

%% Calculando os resíduos de DjPj

for m=1:ord
    for k=1:freq_siz
        mode(m).D(:,:,k)=Ti(:,m,k)*invTi(m,:,k); % matriz Dj idempotentes -> Elimina problemas de indeterminação 
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
opts.passive = 0;
for o = 1:ord
    clear DjPj
    for k = 1:freq_siz
        DjPj(:,:,k) = mode(o).D(:,:,k)*Pj(o,o,k); % Dj.Pj - eq. 3.23 dissertação zanon
        DjPj_t(:,:,k) = transp(DjPj(:,:,k)); % transposta de DjPj
    end
    for k = 1 : numel(Npoles)
        rDjPj = reshape(DjPj_t, size(DjPj_t,2)*size(DjPj_t,1), size(DjPj_t,3)); % reorganiza as dimensões de Yc_t para que fique (Nfases*Nfases) x Namostras % [11; 12; 13; 14; 15; 16; 21; 22; 23; 24; 25; 26; 31; 32; 33; 34; 35; 36 ... 61; 62; 63; 64; 65; 66] vetor coluna
        s = s(:).'; % realizando a transposta de s - 1 x Namostras
          poles = polesPj;
        for j = 1:1
            weight = ones(1,numel(s)); % peso - 1 x Namostras
            [SERH2,poles,~,Happrox,~] = vectfit3(rDjPj, s, poles, weight, opts); % função do VF elaborada por Gustavsen
            Cij2(o) = {SERH2.C};
            if ZYprnt

                linewidth = 3;
                approxColor = ':r';
                for i = 1:ord
                    for j = 1:ord
                        fig_idx=j+fig;
                        data = squeeze(rDjPj(i*j, :));
                        approxData = squeeze(Happrox(i*j, :));
                        figure(fig_idx); clf;
                        for subp = 1:2
                            if subp == 1
                                subplot(2, 1, subp)
                                semilogx(frequency, real(data), 'linewidth', linewidth, 'DisplayName', sprintf('real DjPj_{%d,%d}, mode %d', i,j, o)); hold all;
                                semilogx(frequency, real(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                                legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                            else
                                subplot(2, 1, subp)
                                semilogx(frequency, imag(data), 'linewidth', linewidth, 'DisplayName', sprintf('imag DjPj_{%d,%d}, mode %d ', i,j, o)); hold all;
                                semilogx(frequency, imag(approxData), approxColor, 'linewidth', linewidth, 'DisplayName', sprintf('Fit - %d poles', Npoles)); hold all;
                                legend('-DynamicLegend');xlabel('Frequency [Hz]');ylabel('Magnitude');grid on;

                            end
                        end
                    end
                    fig=fig_idx;
                end
            end
        end
    end
end


%% Escreve os polos e residuos de Yc e H em .dat
%Polos de Yc
polYc = fittedPoles.';
polYc_dupes=find(imag(polYc)<0);
polYc(polYc_dupes)=[];

%Resíduos de Yc
rYc=SERYc.C;

if polYc_dupes~=0         %Exclui os resíduos correspondentes
    for k=1:length(polYc_dupes)
        rYc(:,polYc_dupes(k))=[];
    end
end

%Polos de Pj
polPj2 = polesPj.';
polPj2_dupes=find(imag(polPj2)<0);
polPj2(polPj2_dupes)=[];

m=size(Pj,1);

if polPj2_dupes~=0         %Exclui os resíduos correspondentes
        for i=1:m
            Cij2{1,i}(:,polPj2_dupes)=[];
        end
end

nph=size(Yc,1);
npYc=size(polYc,1);

filename = fullfile(currPath,[jobid '_fitULM.dat']);

fid = fopen(filename,'wt');
fprintf(fid,'%d\n',nph); %numero de fases
fprintf(fid,'%d\n',m); %numero de modos
fprintf(fid,'%d\n',npYc); %numero de polos Yc

for k=1:m
    fprintf(fid,'%d\n',size(polPj2,1)); %numero de polos H para cada modo
end

%taus otimizados
for jj = 1:size(D,1)
    fprintf(fid,'%.16E\n',D(jj));
end

% Escreve os polos de Yc e seus resíduos
for k=1:npYc
    %Escreve o primeiro polo
    pol=polYc(k);
    if isreal(pol)
        if pol>0
            fprintf(fid,'%.16E\n',-(pol)); % Força o polo real a ser negativo
        else
            fprintf(fid,'%.16E\n',pol);
        end
    else
        fprintf(fid,'+%.16E\n%.16E\n',abs(real(pol)),imag(pol)); %Força o polo complexo a ter parte real positiva
    end

    % Escreve os resíduos
    res=rYc(:,k);
    res=reshape(res, [nph nph]);
    res=res.';
    res=triu(res);
    res=reshape(res.', [nph*nph 1]);
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

% Escreve os polos de H e seus resíduos para cada modo
for k=1:m
    for j=1:size(polPj2,1)
         pol=polPj2(j);
        if isreal(pol)
            if pol>0
                fprintf(fid,'%.16E\n',-(pol));
            else
                fprintf(fid,'%.16E\n',pol);
            end
        else
            fprintf(fid,'+%.16E\n%.16E\n',abs(real(pol)),imag(pol));
        end
        
        res = Cij2{1,k}(:,j);
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

% Resíduos de Yc_k0
k0=SERYc.D;
k0=reshape(k0, [nph nph]);
k0=k0.';
k0=triu(k0);
k0=reshape(k0.', [nph*nph 1]);
k0(k0==0)=[];
fprintf(fid,'%.16E\n',real(k0));

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

function [tau_opt, fun] = findoptimtau(f,vel,H,line_length,ERR,opts,Npoles)
j=length(f)-1;
w=2*pi*f;
Ns=length(w);
absH=abs(H);

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
ph=(phase1-phase2);
phase_min=ph; %Phase angle [rad]
tau_mps=(line_length/vel(j)) + phase_min/(w(j));

tmin=1*abs(tau_mps);
tmax=1*abs(line_length/vel(end));
    

options = optimset('Display','none', 'TolX',1e-6);
tau_opt=fminbnd(@(x)fcalc(H,f,x,ERR,opts,Npoles),tmin,tmax,options);
fun = H.*exp(1i.*2.*pi.*f.*tau_opt);

end


function [out] = fcalc(H, f, tau, ERR, opts,Npoles)
current_state=warning;warning('off');
fun = H.*exp(1i.*2.*pi.*f.*tau);
[fit_data,err] = rationalfit(f.',fun.','NPoles',Npoles,'TendsToZero',true,'IterationLimit',[24 40]);
out=err;
warning(current_state);
end
