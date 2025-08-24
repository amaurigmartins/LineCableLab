function [fit_data, ffit, numpol, rmserr] = rationalfit_wrapper_real(fun,freq,opts)

currentWarningState = warning;
warning('off', 'all')

this_var = inputname(1);
rmserr = 1;
j=0;
NpolesMin=opts.Npoles(1);
NpolesMax=opts.Npoles(2);
errtgt = opts.errtgt;
s = 1i*2*pi*freq;
ord = size(fun,2);

while rmserr > errtgt

    numpol=NpolesMin+j;

    [fit_data,errdB] = rationalfit(freq.',fun,'NPoles',numpol,'TendsToZero',opts.TendsToZero,'IterationLimit',opts.Niter);      
    conjugate_par=find(imag(fit_data(1,1).A)~=0);
    if conjugate_par~=0
        for x=1:length(fit_data)
            for k=1:length(conjugate_par)
                par_poles(k,1) = fit_data(1,x).A(conjugate_par(k),1);
                omega_n(k,1) = sqrt(real(par_poles(k,1))^2+imag(par_poles(k,1))^2);
    %             zeta(k,1)= real(par_poles(k,1))/omega_n(k,1);
            end
            for m=1:2:length(conjugate_par)
                proposed_real_(m,1) = real(1i*freq(end)-(omega_n(m,1)-0.01));
                proposed_real_(m+1,1) = real(1i*freq(end)-(omega_n(m+1,1)+0.01));
                fit_data(1,x).A(conjugate_par(m),1) = proposed_real_(m,1);
                fit_data(1,x).A(conjugate_par(m)+1,1) = proposed_real_(m+1,1);
            end
        end
        clear optsR
        optsR.relax = 1;      % use ajuste vetorial com restrição de não trivialidade relaxada
        optsR.stable = 1;     % aplicar pólos estáveis
        optsR.asymp = 2;      % fitting com D~=0, E=0
        optsR.skip_pole = 1;  % pule o cálculo dos polos
        optsR.skip_res = 0;   % não pule o cálculo dos resíduos
        optsR.cmplx_ss = 1;   % crie um modelo de espaço de estado complexo
        optsR.spy1 = 0;       % sem plotagem para o primeiro estágio do VF
        optsR.spy2 = 0;       % criar gráfico de magnitude para ajuste de f(s)

        weight = ones(1,numel(s)); % peso - 1 x Namostras
        [SER,~,rmserr,ffit_,~] = vectfit3(fun.', s,fit_data(1,1).A, weight, optsR);
        for k=1:ord
            fit_data(1,k).C = SER.C(k,:);
            fit_data(1,k).D = 0;
        end
        h_rmserr(j+1,1) = rmserr;
        h_fit_data(j+1,1) = {fit_data};
        for k=1:size(ffit_,1)
                        if opts.plot2
                            figure
                            semilogx(freq,abs(fun(:,k)),'r-','DisplayName','Data','LineWidth',1); hold all
                            semilogx(freq,abs(ffit_(k,:)),'bo','DisplayName',sprintf('Fit, poles = %d',numpol),'LineWidth',1); 
                            semilogx(freq,abs(abs(fun(:,k))-abs(ffit_(k,:).')),'m-','DisplayName','Deviation','LineWidth',1); 
                            legend('-DynamicLegend', 'Location','best')
                            title(sprintf('%s - iteration #%d, rmserr=%1.6f, mode %d',this_var,j,rmserr,k));
                            axis tight
                            grid on;grid minor;
                        end
                    ffit(:,k) = ffit_(k,:);   
        end
        h_ffit(j+1,1) = {ffit};
    else
    
    
    rmserr = 10^(errdB/20);
    h_rmserr(j+1,1) = rmserr;
    h_fit_data(j+1,1) = {fit_data};
    for k=1:length(fit_data)
        [ffit_,~] = freqresp(fit_data(k),freq);
                    if opts.plot2
                        figure
                        semilogx(freq,abs(fun(:,k)),'r-','DisplayName','Data','LineWidth',1); hold all
                        semilogx(freq,abs(ffit_),'bo','DisplayName',sprintf('Fit, poles = %d',numpol),'LineWidth',1); 
                        semilogx(freq,abs(abs(fun(:,k))-abs(ffit_)),'m-','DisplayName','Deviation','LineWidth',1); 
                        legend('-DynamicLegend', 'Location','best')
                        title(sprintf('%s - iteration #%d, rmserr=%1.6f, mode %d',this_var,j,rmserr,k));
                        axis tight
                        grid on;grid minor;
                    end
         ffit(:,k) = ffit_;   
    end
    h_ffit(j+1,1) = {ffit};
    
    end
    
    if rmserr < errtgt
        break;
    end
    
    if opts.err_par
        if j > 0
            noErrorImprovement = (abs( h_rmserr(end)-h_rmserr(end-1)) < errtgt ); %adding poles does not improve the response
            addingPolIncreasesError = ( h_rmserr(end) > h_rmserr(end-1) ); %adding poles makes things worse >:|
            if noErrorImprovement || addingPolIncreasesError
                fit_data=h_fit_data{end-1};
                numpol=length(fit_data(1,1).A);
                rmserr = h_rmserr(end-1);
                ffit = h_ffit{end-1};
                break
            end
        end
    end
    if numpol == NpolesMax
        break;
    end         % Reached maximum number of poles
    j=j+1;
end
  warning(currentWarningState)

end
