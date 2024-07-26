function [fit_data, ffit, numpol, rmserr] = rationalfit_wrapper_real(fun,freq,opts)

currentWarningState = warning;
warning('off', 'all')

this_var = inputname(1);
rmserr = 1;
i=0;
NpolesMin=opts.Npoles(1);
NpolesMax=opts.Npoles(2);
errtgt = opts.errtgt;
tolThreshold =errtgt;
Ns=length(freq);
s = j*2*pi*freq;

while rmserr > errtgt

    numpol=NpolesMin+i;

    [fit_data,errdB] = rationalfit(freq.',fun,'NPoles',numpol,'TendsToZero',opts.TendsToZero,'IterationLimit',opts.Niter);      
    if i==0
        if ~isreal(fit_data(1,1).A)
              warning('Rationalfit was unable to fit the data with real poles only. ATP results might be unstable or inaccurate.');
              warning('Will try to fit Zthe data using Bode method. Beware possible innaccuracies.');
              tol = opts.tol; %Tolerance in decibels to set a new pole and/or new zero
              try
                    [P,Z,k,errdB] = Bode_process(abs(fun),freq,Ns,tol);
                    as = k.*poly(Z); bs = poly(P); % Polynomials
                    [r,p,ks] = residue(as,bs); % Poles, residues and constant term
              catch
                  error('Bode method failed to fit one or more modes. No PCH file will be written.');
              end
              TF=isempty(ks);if(TF==1);ks=0;end

              fit_data.A = p;
              fit_data.C = r;
              fit_data.D = ks;
              numpol = length(p);
              rmserr = 10^(errdB/20);
              h_rmserr(i+1,1) = rmserr;
              h_fit_data(i+1,1) = {fit_data};
              for k=1:length(fit_data)
                    [ffit_,~] = freqresp(fit_data(k),freq);
                    freq_extrap=logspace(ceil(log10(freq(end))),ceil(log10(freq(end)))+4,80).';
                    freq_extrap = [freq(1:end-1); freq_extrap];
                    [ffit_extrap,~] = freqresp(fit_data(k),freq_extrap);
                        if opts.plot2
                            figure
                            semilogx(freq,abs(fun(:,k)),'r-','DisplayName','Data','LineWidth',1); hold all
                            semilogx(freq_extrap,abs(ffit_extrap),'bo','DisplayName',sprintf('Fit, poles = %d',numpol),'LineWidth',1); 
                            semilogx(freq,abs(abs(fun(:,k))-abs(ffit_)),'m-','DisplayName','Deviation','LineWidth',1); 
                            legend('-DynamicLegend', 'Location','best')
                            title(sprintf('%s - iteration #%d, rmserr=%1.6f, mode %d',this_var,i,rmserr,k));
                            axis tight
                            grid on;grid minor;
                        end
                    ffit(:,k) = ffit_;   
              end
              h_ffit(i+1,1) = {ffit};
              break
        end
    end
    
    rmserr = 10^(errdB/20);
    h_rmserr(i+1,1) = rmserr;
    h_fit_data(i+1,1) = {fit_data};
    for k=1:length(fit_data)
                [ffit_,~] = freqresp(fit_data(k),freq);
                freq_extrap=logspace(ceil(log10(freq(end))),ceil(log10(freq(end)))+4,80).';
                freq_extrap = [freq(1:end-1); freq_extrap];
                [ffit_extrap,~] = freqresp(fit_data(k),freq_extrap);
                    if opts.plot2
                        figure
                        semilogx(freq,abs(fun(:,k)),'r-','DisplayName','Data','LineWidth',1); hold all
                        semilogx(freq_extrap,abs(ffit_extrap),'bo','DisplayName',sprintf('Fit, poles = %d',numpol),'LineWidth',1); 
                        semilogx(freq,abs(abs(fun(:,k))-abs(ffit)),'m-','DisplayName','Deviation','LineWidth',1); 
                        legend('-DynamicLegend', 'Location','best')
                        title(sprintf('%s - iteration #%d, rmserr=%1.6f, mode %d',this_var,i,rmserr,k));
                        axis tight
                        grid on;grid minor;
                    end
                ffit(:,k) = ffit_;   
    end
    h_ffit(i+1,1) = {ffit};
    

    if rmserr < errtgt
        break;
    end
    if i > 0
        noErrorImprovement = (abs( h_rmserr(end)-h_rmserr(end-1)) < tolThreshold ); %adding poles does not improve the response
        addingPolIncreasesError = ( h_rmserr(end) > h_rmserr(end-1) ); %adding poles makes things worse >:|
        if noErrorImprovement || addingPolIncreasesError
            fit_data=h_fit_data{end-1};
            numpol=length(fit_data(1,1).A);
            rmserr = h_rmserr(end-1);
            ffit = h_ffit{end-1};
            break
        end
        if ~isreal(fit_data(1,1).A)
            fit_data=h_fit_data{end-1};
            numpol=length(fit_data(1,1).A);
            rmserr = h_rmserr(end-1);
            ffit = h_ffit{end-1};
            break;
        end
    end

    if numpol == NpolesMax
        break;
    end         % Reached maximum number of poles
    i=i+1;
end
  warning(currentWarningState)

end
