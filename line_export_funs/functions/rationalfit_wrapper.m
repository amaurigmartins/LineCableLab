function [fit_data, ffit, numpol, rmserr] = rationalfit_wrapper(fun,freq,opts)

currentWarningState = warning;
warning('off', 'all')

this_var = inputname(1);
rmserr = 1;
i=0;
NpolesMin=opts.Npoles(1);
NpolesMax=opts.Npoles(2);
errtgt = opts.errtgt;
tolThreshold =errtgt;

while rmserr > errtgt

    numpol=NpolesMin+i;

    [fit_data,errdB] = rationalfit(freq.',fun,'NPoles',numpol,'TendsToZero',opts.TendsToZero,'IterationLimit',opts.Niter);
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
                semilogx(freq,abs(fun(:,k)),'r-','DisplayName','data','LineWidth',1); hold all
                semilogx(freq_extrap,abs(ffit_extrap),'bo','DisplayName',sprintf('fit, poles = %d',numpol),'LineWidth',1); 
                semilogx(freq,abs(fun(:,k)-ffit_),'m-','DisplayName','deviation','LineWidth',1); 
                xlim([min(freq) max(freq_extrap)]);
                legend('-DynamicLegend', 'Location','best')
                title(sprintf('%s - iteration #%d, rmserr=%1.6f, mode %d',this_var,i,rmserr,k))
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
    end

    if numpol == NpolesMax
        break;
    end         % Reached maximum number of poles
    i=i+1;
end %of main loop

warning(currentWarningState)

% backup
% currentWarningState = warning;
% warning('off', 'all')
% 
% this_var = inputname(1);
% rmserr = 1;
% i=0;
% NpolesMin=opts.Npoles(1);
% NpolesMax=opts.Npoles(2);
% errtgt = opts.errtgt;
% tolThreshold =errtgt;
% 
% while rmserr > errtgt
% 
%     numpol=NpolesMin+i;
% 
%     [fit_data,errdB] = rationalfit(freq.',fun,'NPoles',numpol,'TendsToZero',opts.TendsToZero,'IterationLimit',opts.Niter);
%     rmserr = 10^(errdB/20);
%     h_fit_data(i+1,1) = {fit_data};
%     [ffit,~] = freqresp(fit_data,freq);
%     freq_extrap=logspace(ceil(log10(freq(end))),ceil(log10(freq(end)))+4,80).';
%     freq_extrap = [freq(1:end-1); freq_extrap];
%     [ffit_extrap,~] = freqresp(fit_data,freq_extrap);
%     h_rmserr(i+1,1) = rmserr;
% 
%     if opts.plot2
%         figure
%         % real_poles=imag(fit_data.A)==0;
%         % [zero_freq, pol_freq, ~, ~] = zpk(fit_data);
%         % pol_val = interp1(freq_extrap, abs(ffit_extrap), pol_freq);
%         % zero_val = interp1(freq_extrap, abs(ffit_extrap), zero_freq);
%         semilogx(freq,abs(fun),'r-','DisplayName','data','LineWidth',1); hold all
%         semilogx(freq_extrap,abs(ffit_extrap),'bo','DisplayName',sprintf('fit, poles = %d',numpol),'LineWidth',1); 
%         semilogx(freq,abs(fun-ffit.'),'m-','DisplayName','deviation','LineWidth',1); 
%         % semilogx(pol_freq(real_poles),pol_val(real_poles),'^','MarkerSize',8, 'MarkerFaceColor','blue','DisplayName','real pole')
%         % semilogx(pol_freq(~real_poles),pol_val(~real_poles),'*','MarkerSize',8, 'MarkerFaceColor','red','DisplayName','complex pole')
%         % semilogx(zero_freq,zero_val,'o','MarkerSize',8, 'MarkerFaceColor','green','DisplayName','zero')
%         xlim([min(freq) max(freq_extrap)]);
%         legend('-DynamicLegend', 'Location','best')
%         title(sprintf('%s fit - iteration #%d, rmserr=%1.6f',this_var,i,rmserr))
%     end
% 
%     if rmserr < errtgt
%         break;
%     end
%     if i > 0
%         noErrorImprovement = (abs( h_rmserr(end)-h_rmserr(end-1)) < tolThreshold ); %adding poles does not improve the response
%         addingPolIncreasesError = ( h_rmserr(end) > h_rmserr(end-1) ); %adding poles makes things worse >:|
%         if noErrorImprovement || addingPolIncreasesError
%             % warning('Rationalfit stopped by "No Error Improvement" or "Adding Poles Increases Error"')
%             fit_data=h_fit_data{end-1};
%             numpol=length(fit_data.A);
%             [ffit,~] = freqresp(fit_data,freq);
%             rmserr = h_rmserr(end-1);
%             break
%         end
%     end
% 
%     if numpol == NpolesMax
%         break;
%     end         % Reached maximum number of poles
%     i=i+1;
% end
% 
% warning(currentWarningState)

end %of function