function [pol, res, infval, N, fit, rmserr] = vectfit_wrapper(f,freq,errtgt,opts)

if size(f,1) > 1
    f=f.';
end

if size(freq,1) > 1
    freq=freq.';
end

Ns=length(freq);

warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:rankDeficientMatrix')

s=1i*2*pi*freq;
if opts.weightscheme==1
    weight=ones(1,Ns);
elseif opts.weightscheme==2
    weight=zeros(1,Ns);
    for k=1:Ns
        weight(1,k)=1/sqrt(norm(freq(:,k)));
    end
elseif opts.weightscheme==3
    weight=zeros(1,Ns);
    for k=1:Ns
        weight(1,k)=1/(norm(freq(:,k)));
    end
end

%Initial poles for Vector Fitting:
if isfield(opts,'NORD')
    N=opts.NORD;
else
    N=1; %order of approximation
end

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
% in=log(imag(s(1))/(2*pi))/log(10);
% fin=log(imag(s(end))/(2*pi))/log(10);
% poles=-2*pi*logspace(in,fin,N); %Initial poles

if isfield(opts,'Niter')
    for k=1:opts.Niter
        [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
    end    
else
    MAXORDER=40;
    TOL = 1e-6;
    MAXITER = 100;
    rmserr = 1;
    i=0;
    if opts.output_messages
        fprintf('Starting fitting...\n')
    end
    while rmserr > errtgt
        [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
        [R,a]=ss2pr(SER.A,SER.B,SER.C);
        SER.R=R;
        SER.poles=a;
        if opts.passive==1
            % Passivity enforcement
            % clear opts;
            opts2.Niter_out = 100;
            opts2.Niter_in = 100;
            opts2.parametertype='Y';
            opts2.cmplx_ss=1;
            opts2.outputlevel = 1;
            [SER,~,~]=RPdriver(SER,s,opts2);
        end


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
res=squeeze(R);
pol=a;
infval = SER.D;
N = length(pol);

warning('on', 'MATLAB:nearlySingularMatrix')
warning('on', 'MATLAB:rankDeficientMatrix')

end