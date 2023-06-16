function [] = plotwavchar_fun(freq,ord,a,vel,Zch_mod,Ti,sigma_g_total,erg_total,tag)

%% Plot attenuation, velocity, Zch_mod
e0=8.85418782*1e-12;
cpcond= sigma_g_total + (1i .* 2 .* pi .* freq .* erg_total .* e0);
% vel=abs(vel);

% set_plot_params()
figure('Name', ['AtnConst_PhaseVel_' tag])

subplot(2,1,1)
for o=1:ord
    loglog(freq,a(:,o),'LineWidth',2);
    hold all
    lgd{o} = strcat('mode # ',num2str(o));
end
xlabel('Frequency [Hz]');
ylabel('\alpha [Np/m]');
title('Attenuation constant')
legend(lgd);
grid on;

hold off
subplot(2,1,2)
for o=1:ord
    semilogx(freq,vel(:,o),'LineWidth',2)
    hold all
    lgd{o} = strcat('mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('\upsilon [m/s]')
title('Phase velocity')
legend(lgd);
grid on

if ~all(sigma_g_total == sigma_g_total(1)) % only plot this if sigma is not constant
    set_plot_params()
    figure('Name', ['AtnConst_PhaseVel_EarthResis_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(1./sigma_g_total,a(:,o),'LineWidth',2);
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('\alpha [Np/m]');
    title('Attenuation constant')
    axis tight
    legend(lgd);
    grid on;
    
    hold off
    subplot(2,1,2)
    for o=1:ord
        plot(1./sigma_g_total,vel(:,o),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('\upsilon [m/s]')
    title('Phase velocity')
    axis tight
    legend(lgd);
    grid on
    
    set_plot_params()
    figure('Name', ['AtnConst_PhaseVel_EarthPerm_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(erg_total,a(:,o),'LineWidth',2);
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    xlabel('Relative earth permittivity [p.u.]')
    ylabel('\alpha [Np/m]');
    title('Attenuation constant')
    axis tight
    legend(lgd);
    grid on;
    
    hold off
    subplot(2,1,2)
    for o=1:ord
        plot(erg_total,vel(:,o),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('Relative earth permittivity [p.u.]')
    ylabel('\upsilon [m/s]')
    title('Phase velocity')
    axis tight
    legend(lgd);
    grid on
    
    
   
    
    
    set_plot_params()
    figure('Name', ['AtnConst_PhaseVel_ComplexEarthResis_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(abs(1./cpcond),a(:,o),'LineWidth',2);
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    xlabel('abs(1/(\sigma + j \omega \epsilon)) [\Omega.m]')
    ylabel('\alpha [Np/m]');
    title('Attenuation constant')
    axis tight
    legend(lgd);
    grid on;
    
    hold off
    subplot(2,1,2)
    for o=1:ord
        plot(abs(1./cpcond),vel(:,o),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('abs(1/(\sigma + j \omega \epsilon)) [\Omega.m]')
    ylabel('\upsilon [m/s]')
    title('Phase velocity')
    axis tight
    legend(lgd);
    grid on
    
end

set_plot_params()
figure('Name', ['CharImped_' tag])

subplot(2,1,1)
for o=1:ord
    semilogx(freq,abs(Zch_mod(:,o)),'LineWidth',2)
    hold all
    lgd{o} = strcat('mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('Magnitude [\Omega]')
title('Characteristic impedance')
legend(lgd);
grid on

subplot(2,1,2)
for o=1:ord
    semilogx(freq,radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
    hold all
    lgd{o} = strcat('mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('Angle [deg]')
legend(lgd);
grid on

if ~all(sigma_g_total == sigma_g_total(1)) % only plot this if sigma is not constant
    set_plot_params()
    figure('Name', ['CharImpedEarthResis_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(1./sigma_g_total,abs(Zch_mod(:,o)),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('Magnitude [\Omega]')
    title('Characteristic impedance')
    axis tight
    legend(lgd);
    grid on
    
    subplot(2,1,2)
    for o=1:ord
        plot(1./sigma_g_total,radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('Angle [deg]')
    axis tight
    legend(lgd);
    grid on
    
    set_plot_params()
    figure('Name', ['CharImpedEarthPerm_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(erg_total,abs(Zch_mod(:,o)),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('Relative earth permittivity [p.u.]')
    ylabel('Magnitude [\Omega]')
    title('Characteristic impedance')
    axis tight
    legend(lgd);
    grid on
    
    subplot(2,1,2)
    for o=1:ord
        plot(erg_total,radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('Relative earth permittivity [p.u.]')
    ylabel('Angle [deg]')
    axis tight
    legend(lgd);
    grid on
    
 
    
    set_plot_params()
    figure('Name', ['CharImpedComplexEarthResis_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(abs(1./cpcond),abs(Zch_mod(:,o)),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('abs(1/(\sigma + j \omega \epsilon)) [\Omega.m]')
    ylabel('Magnitude [\Omega]')
    title('Characteristic impedance')
    axis tight
    legend(lgd);
    grid on
    
    subplot(2,1,2)
    for o=1:ord
        plot(abs(1./cpcond),radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
        hold all
        lgd{o} = strcat('mode # ',num2str(o));
    end
    hold off
    xlabel('abs(1/(\sigma + j \omega \epsilon)) [\Omega.m]')
    ylabel('Angle [deg]')
    axis tight
    legend(lgd);
    grid on
    
end

%[sk col k_scaling w_corr h_corr k_width_height fnt_scaling
set_plot_params([1 2 2 .8 .6 1 .7])
figure('Name', ['TransfMat_' tag])
%set_plot_params(.75,2.5,1.15) % a lot of trial n error to get this almost right
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,abs(Ti(:,k+(o-1)*ord)),'LineWidth',2)
        %semilogx(freq,abs(Ti(:,k+(o-1)*ord)),freq,abs(Ti_vf(:,k+(o-1)*ord))) % ??? vector fitting
        hold all
        if o==1;lgd{o} = strcat('mode # ',num2str(o));end
    end
    legend(lgd);
    grid on
    xlabel('Frequency [Hz]')
    if k==1;title('Transformation matrix magnitude');end
    hold off
    
    subplot(ord,2,2*k)
    for o=1:ord
        %semilogx(freq,radtodeg(angle(Ti(:,k+(o-1)*ord))))
        semilogx(freq,radtodeg(unwrap(angle(Ti(:,k+(o-1)*ord)))),'LineWidth',2)
        if o==1;lgd{o} = strcat('mode # ',num2str(o));end
        %semilogx(freq,radtodeg(acos(real(Ti(:,k+(o-1)*ord))./abs(Ti(:,k+(o-1)*ord)))))
        %semilogx(freq,radtodeg(asin(imag(Ti(:,k+(o-1)*ord))./abs(Ti(:,k+(o-1)*ord)))))
        %semilogx(freq,radtodeg(unwrap(angle(Ti(:,k+(o-1)*ord)))),freq,radtodeg(unwrap(angle(Ti_vf(:,k+(o-1)*ord))))) % ??? vector fitting
        hold all
    end
%     legend(lgd);
    grid on
    xlabel('Frequency [Hz]')
    if k==1;title('Transformation matrix angle');end
    hold off
end