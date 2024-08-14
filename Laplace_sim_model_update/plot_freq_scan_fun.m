function [] = plot_freq_scan_fun(V,ord,freq,tag)

    Vmag=abs(V);

    for i=1:ord
        plottitle{i}=sprintf('Phase #%d - Sending terminal',i);
        plottitle{i+ord}=sprintf('Phase #%d - Receiving terminal',i);
    end

    for o=1:2*ord
        figure('Name', ['FreqScan' num2str(o) '_' tag])
        semilogx(freq,Vmag(o,:),'LineWidth',2);
        xlabel('Frequency [Hz]')
        ylabel('Magnitude [pu]')
        title(plottitle{o});
        grid on
    end

end