function [pol, res, infval, N, fit, rmserr] = vectfit_wrapper(f,freq,errtgt,opts)

if size(f,1) > 1
    f=f';
end

if size(freq,1) > 1
    freq=freq';
end

Ns=length(freq);

opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %Do NOT skip identification of residues (C,D,E)
opts.cmplx_ss=1;   %Create real state space model?


opts.spy1=0;       %No plotting for first stage of vector fitting
opts.logx=1;       %Use logarithmic abscissa axis
opts.logy=0;       %Use logarithmic ordinate axis
opts.errplot=0;    %Include deviation in magnitude plot
opts.phaseplot=0;  %Also produce plot of phase angle (in addition to magnitiude)
opts.legend=0;     %Do include legends in plots


addpath('vfit3')

warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:rankDeficientMatrix')

s=1i*2*pi*freq;
if opts.weightscheme==1
    weight=ones(1,Ns);
elseif opts.weightscheme==2
    weight=zeros(1,Ns);
    for k=1:Ns
        weight(1,k)=1/sqrt(norm(f(:,k)));
    end
end

%Initial poles for Vector Fitting:
N=1; %order of approximation

if opts.firstguesstype == 1
    poles=-2*pi*logspace(log10(freq(1)),log10(freq(end)),N); %Initial poles
elseif opts.firstguesstype == 2
    % Complex conjugate pairs, logarithmically spaced :
    bet=logspace(log10(freq(1)),log10(freq(end)),N/2);
    poles=[];
    for n=1:length(bet)
        alf=-bet(n)*1e-2;
        poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ];
    end
elseif opts.firstguesstype == 3
    % Complex conjugate pairs, mix. linearly and logarithmically spaced :
    nu=0.001;
    bet=linspace(s(1)/1i,s(Ns)/1i,ceil((N-1)/4));
    poles1=[];
    for n=1:length(bet)
        alf=-nu*bet(n);
        poles1=[poles1 (alf-1i*bet(n)) (alf+1i*bet(n)) ];
    end
    bet=logspace(log10(s(1)/1i),log10(s(Ns)/1i),2+floor(N/4));
    bet(1)=[];bet(end)=[];
    poles2=[];
    for n=1:length(bet)
        alf=-nu*bet(n);
        poles2=[poles2 (alf-1i*bet(n)) (alf+1i*bet(n)) ];
    end
    poles=[poles1 poles2];
end

MAXORDER=60;
TOL = 1e-6;
MAXITER = 100;
rmserr = 1;
i=0;
if opts.output_messages
    fprintf('Starting fitting with complex poles...\n')
end
while rmserr > errtgt
    [SER,poles,rmserr,fit]=vectfit4(f,s,poles,weight,opts);
    rmserr = abs(rmserr / mean(f));
    i=i+1;
    if opts.output_messages
        fprintf('Iteration #%d...\n', i)
        fprintf('--- Relative RMS error is %1.4f%%\n', rmserr*100)
    end
    err_hist(i)=rmserr;
    if i>=2
        dE = err_hist(i)-err_hist(i-1);
        if dE < TOL
            N = N+1;
            if opts.firstguesstype == 1
                poles=-2*pi*logspace(log10(freq(1)),log10(freq(end)),N); %Initial poles
            elseif opts.firstguesstype == 2
                % Complex conjugate pairs, logarithmically spaced :
                bet=logspace(log10(freq(1)),log10(freq(end)),N/2);
                poles=[];
                for n=1:length(bet)
                    alf=-bet(n)*1e-2;
                    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ];
                end
            elseif opts.firstguesstype == 3
                % Complex conjugate pairs, mix. linearly and logarithmically spaced :
                nu=0.001;
                bet=linspace(s(1)/1i,s(Ns)/1i,ceil((N-1)/4));
                poles1=[];
                for n=1:length(bet)
                    alf=-nu*bet(n);
                    poles1=[poles1 (alf-1i*bet(n)) (alf+1i*bet(n)) ];
                end
                bet=logspace(log10(s(1)/1i),log10(s(Ns)/1i),2+floor(N/4));
                bet(1)=[];bet(end)=[];
                poles2=[];
                for n=1:length(bet)
                    alf=-nu*bet(n);
                    poles2=[poles2 (alf-1i*bet(n)) (alf+1i*bet(n)) ];
                end
                poles=[poles1 poles2];
            end
            if opts.output_messages
                if rmserr > errtgt
                    fprintf('Error diff reached TOL without meeting errtgt. Increasing fit order.\n')
                end
            end
        end
    end
    if i == MAXITER
        if opts.output_messages
            fprintf('MAXITER reached. Stopping now.\n')
        end
        break
    end
    if N == MAXORDER
        if opts.output_messages
            fprintf('MAXORDER reached. Stopping now.\n')
        end
        break
    end
end


if opts.output_messages
    fprintf('\n********************************************************************\n')
    fprintf('Process converged after iteration #%d...\n', i)
    fprintf('--- Relative RMS error is %1.4f%%\n', rmserr*100)
    fprintf('--- Fit order is %d\n', N)
end

[R,a]=ss2pr(SER.A,SER.B,SER.C);
SER.R=R;
SER.poles=a;

if opts.passive==1
    % Passivity enforcement
    % clear opts;
    opts.Niter_out = 100;
    opts.Niter_in = 100;
    opts.parametertype='Y';
    opts.cmplx_ss=1;
    opts.outputlevel = 1;
    [SER,~,~]=RPdriver(SER,s,opts);
end
% SER = reducecmplx(SER);

% [R,a]=ss2pr(SER.A,SER.B,SER.C);
% % Remove any zero residues, if any
% R=squeeze(R);
% a(R==0)=[];
% R(R==0)=[];
% SER.R=R;
% SER.poles=a;

% %Final stage to get actual fitted values
[SER,poles,rmserr,fit]=vectfit4(f,s,poles,weight,opts);

[R,a]=ss2pr(SER.A,SER.B,SER.C);
res=R;
pol=a;
infval = SER.D;
N = length(pol);

warning('on', 'MATLAB:nearlySingularMatrix')
warning('on', 'MATLAB:rankDeficientMatrix')

end