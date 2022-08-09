function [] = plotwavchar_fun(freq,ord,g,Zch_mod,Ti)

%% Plot attenuation, velocity, Zch_mod
a=real(g);
b=imag(g);
freq_siz=length(freq);
vel=zeros(freq_siz,ord);
for o=1:ord
    vel(:,o)=(2*pi*freq)./b(:,o);
end

figure()
subplot(2,1,1)
for o=1:ord
    loglog(freq,a(:,o),'LineWidth',2);
    hold all
end
xlabel('frequency (Hz)');
ylabel('\alpha (Np/m)');
title('attenuation constant');
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
xlabel('frequency (Hz)')
ylabel('\upsilon (m/s)')
title('velocity')
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

figure()
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,abs(Ti(:,k+(o-1)*ord)),'LineWidth',2)
        %semilogx(freq,abs(Ti(:,k+(o-1)*ord)),freq,abs(Ti_vf(:,k+(o-1)*ord))) % ??? vector fitting
        grid;
        hold all
    end
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
    xlabel('frequency (Hz)')
    title('Transformation matrix')
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
    xlabel('frequency (Hz)')
    title('Transformation matrix')
    hold off
end



figure()
subplot(2,1,1)
for o=1:ord
    semilogx(freq,abs(Zch_mod(:,o)),'LineWidth',2)
    hold all
end
hold off
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega)')
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
xlabel('frequency (Hz)')
ylabel('Angle (degrees)')
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
