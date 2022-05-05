function [] = plotwavchar_fun(freq,ord,a,vel,Zch_mod,Ti,sigma_g_total,erg_total,tag)

%% Plot attenuation, velocity, Zch_mod

set_plot_params()
figure('Name', ['AtnConst_PhaseVel_' tag])

subplot(2,1,1)
for o=1:ord
    loglog(freq,a(:,o),'LineWidth',2);
    hold all
end
xlabel('Frequency [Hz]');
ylabel('\alpha [Np/m]');
title('Attenuation constant')
if ord==3
    legend('mode #1','mode #2','mode #3')
end
if ord==4
    legend('mode #1','mode #2','mode #3','mode #4')
end
if ord==5
    legend('mode #1','mode #2','mode #3','mode #4','mode #5')
end
if ord==6
    legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
end
grid;

hold off
subplot(2,1,2)
for o=1:ord
    semilogx(freq,vel(:,o),'LineWidth',2)
    hold all
end
hold off
xlabel('Frequency [Hz]')
ylabel('\upsilon [m/s]')
title('Phase velocity')

if ord==3
    legend('mode #1','mode #2','mode #3')
end
if ord==4
    legend('mode #1','mode #2','mode #3','mode #4')
end
if ord==5
    legend('mode #1','mode #2','mode #3','mode #4','mode #5')
end
if ord==6
    legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
end
grid

if ~all(sigma_g_total == sigma_g_total(1)) % only plot this if sigma is not constant
    set_plot_params()
    figure('Name', ['AtnConst_PhaseVel_EarthResis_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(1./sigma_g_total,a(:,o),'LineWidth',2);
        hold all
    end
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('\alpha [Np/m]');
    title('Attenuation constant')
    axis tight
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid;
    
    hold off
    subplot(2,1,2)
    for o=1:ord
        plot(1./sigma_g_total,vel(:,o),'LineWidth',2)
        hold all
    end
    hold off
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('\upsilon [m/s]')
    title('Phase velocity')
    axis tight
    
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid
    
    set_plot_params()
    figure('Name', ['AtnConst_PhaseVel_EarthPerm_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(erg_total,a(:,o),'LineWidth',2);
        hold all
    end
    xlabel('Earth permittivity [F/m]')
    ylabel('\alpha [Np/m]');
    title('Attenuation constant')
    axis tight
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid;
    
    hold off
    subplot(2,1,2)
    for o=1:ord
        plot(erg_total,vel(:,o),'LineWidth',2)
        hold all
    end
    hold off
    xlabel('Earth permittivity [F/m]')
    ylabel('\upsilon [m/s]')
    title('Phase velocity')
    axis tight
    
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid
end

set_plot_params()
figure('Name', ['CharImped_' tag])

subplot(2,1,1)
for o=1:ord
    semilogx(freq,abs(Zch_mod(:,o)),'LineWidth',2)
    hold all
end
hold off
xlabel('Frequency [Hz]')
ylabel('Magnitude [\Omega]')
title('Characteristic impedance')
if ord==3
    legend('mode #1','mode #2','mode #3')
end
if ord==4
    legend('mode #1','mode #2','mode #3','mode #4')
end
if ord==5
    legend('mode #1','mode #2','mode #3','mode #4','mode #5')
end
if ord==6
    legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
end
grid

subplot(2,1,2)
for o=1:ord
    semilogx(freq,radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
    hold all
end
hold off
xlabel('Frequency [Hz]')
ylabel('Angle [deg]')
if ord==3
    legend('mode #1','mode #2','mode #3')
end
if ord==4
    legend('mode #1','mode #2','mode #3','mode #4')
end
if ord==5
    legend('mode #1','mode #2','mode #3','mode #4','mode #5')
end
if ord==6
    legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
end
grid

if ~all(sigma_g_total == sigma_g_total(1)) % only plot this if sigma is not constant
    set_plot_params()
    figure('Name', ['CharImpedEarthResis_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(1./sigma_g_total,abs(Zch_mod(:,o)),'LineWidth',2)
        hold all
    end
    hold off
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('Magnitude [\Omega]')
    title('Characteristic impedance')
    axis tight
    
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid
    
    subplot(2,1,2)
    for o=1:ord
        plot(1./sigma_g_total,radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
        hold all
    end
    hold off
    xlabel('Earth resistivity [\Omega.m]')
    ylabel('Angle [deg]')
    axis tight
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid
    
    set_plot_params()
    figure('Name', ['CharImpedEarthPerm_' tag])
    
    subplot(2,1,1)
    for o=1:ord
        plot(erg_total,abs(Zch_mod(:,o)),'LineWidth',2)
        hold all
    end
    hold off
    xlabel('Earth permittivity [F/m]')
    ylabel('Magnitude [\Omega]')
    title('Characteristic impedance')
    axis tight
    
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid
    
    subplot(2,1,2)
    for o=1:ord
        plot(erg_total,radtodeg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
        hold all
    end
    hold off
    xlabel('Earth permittivity [F/m]')
    ylabel('Angle [deg]')
    axis tight
    if ord==3
        legend('mode #1','mode #2','mode #3')
    end
    if ord==4
        legend('mode #1','mode #2','mode #3','mode #4')
    end
    if ord==5
        legend('mode #1','mode #2','mode #3','mode #4','mode #5')
    end
    if ord==6
        legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
    end
    grid
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
        grid;
        hold all
    end
    if k==1
        if ord==3
            legend('mode #1','mode #2','mode #3')
        end
        if ord==4
            legend('mode #1','mode #2','mode #3','mode #4')
        end
        if ord==5
            legend('mode #1','mode #2','mode #3','mode #4','mode #5')
        end
        if ord==6
            legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
        end
    end
    xlabel('Frequency [Hz]')
    if k==1;title('Transformation matrix magnitude');end
    hold off
    
    subplot(ord,2,2*k)
    for o=1:ord
        %semilogx(freq,radtodeg(angle(Ti(:,k+(o-1)*ord))))
        semilogx(freq,radtodeg(unwrap(angle(Ti(:,k+(o-1)*ord)))),'LineWidth',2)
        grid;
        %semilogx(freq,radtodeg(acos(real(Ti(:,k+(o-1)*ord))./abs(Ti(:,k+(o-1)*ord)))))
        %semilogx(freq,radtodeg(asin(imag(Ti(:,k+(o-1)*ord))./abs(Ti(:,k+(o-1)*ord)))))
        %semilogx(freq,radtodeg(unwrap(angle(Ti(:,k+(o-1)*ord)))),freq,radtodeg(unwrap(angle(Ti_vf(:,k+(o-1)*ord))))) % ??? vector fitting
        hold all
    end
%     if k==1
%         if ord==3
%             legend('mode #1','mode #2','mode #3')
%         end
%         if ord==4
%             legend('mode #1','mode #2','mode #3','mode #4')
%         end
%         if ord==5
%             legend('mode #1','mode #2','mode #3','mode #4','mode #5')
%         end
%         if ord==6
%             legend('mode #1','mode #2','mode #3','mode #4','mode #5','mode #6')
%         end
%     end
    xlabel('Frequency [Hz]')
    if k==1;title('Transformation matrix angle');end
    hold off
end