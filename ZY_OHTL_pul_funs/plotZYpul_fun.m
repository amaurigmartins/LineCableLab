function [] = plotZYpul_fun(f,ord,data,jobid,opts)

if nargin==4
    opts.cartesian=false;
    opts.loglog=false;
    opts.plotall=false;
else
    if ~isfield(opts,'cartesian'); opts.cartesian=false; end
    if ~isfield(opts,'loglog'); opts.loglog=false; end
    if ~isfield(opts,'plotall'); opts.plotall=false; end
end

numvars=length(data);
if opts.plotall
    irange=ord;
else
    irange=1;
end

if startsWith(data(1).VarName,'Ztot')
    varSymbol='Z';
    varQty='impedance';
    varUnit='[\Omega/m]';
elseif startsWith(data(1).VarName,'Ytot')
    varSymbol='Y';
    varQty='admittance';
    varUnit='[Sm]';
end

% Plot all results
for i=1:irange
    for j=1:ord
        if j==i; plotLabel='Self';else;plotLabel='Mutual';end
        figure('Name', sprintf('%s%s%d%d_%s',plotLabel,varSymbol,i,j,jobid))
        for k=1:numvars
            subplot(2,1,1)
            if opts.cartesian
                val=squeeze(real(data(k).Values(i,j,:)));
                ylbl=['Real part ' varUnit];
            else
                val=squeeze(abs(data(k).Values(i,j,:)));
                ylbl=['Magnitude ' varUnit];
            end
            lgdLabel=data(k).Label;
            if opts.loglog
                loglog(f,val,'LineWidth',2,'DisplayName',lgdLabel);hold all
            else
                semilogx(f,val,'LineWidth',2,'DisplayName',lgdLabel);hold all
            end
        end
        legend('-DynamicLegend','Location','best');
        ylabel(sprintf('%s',ylbl))
        xlabel('Frequency [Hz]')
        grid on
        title(sprintf('%s %s - %s%d%d',plotLabel,varQty,varSymbol,i,j))

        for k=1:numvars
            subplot(2,1,2)

            if opts.cartesian
                val=squeeze(imag(data(k).Values(i,j,:)));
                ylbl=['Imaginary part ' varUnit];
                if opts.loglog
                    loglog(f,val,'LineWidth',2,'DisplayName',lgdLabel);hold all
                else
                    semilogx(f,val,'LineWidth',2,'DisplayName',lgdLabel);hold all
                end
            else
                val=squeeze(rad2deg(unwrap(squeeze(angle(data(k).Values(i,j,:))))));
                ylbl='Angle [º]';
                semilogx(f,val,'LineWidth',2,'DisplayName',lgdLabel);hold all
            end
            lgdLabel=data(k).Label;

        end
        ylabel(sprintf('%s',ylbl))
        xlabel('Frequency [Hz]')
        grid on
    end
end

end