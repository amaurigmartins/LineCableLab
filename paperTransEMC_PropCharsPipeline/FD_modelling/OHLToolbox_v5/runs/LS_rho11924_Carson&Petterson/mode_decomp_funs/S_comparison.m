%S_Ratio=S_Theofilos./S_Sunde;
S_Ratio=S_Sunde./S_Theofilos;

figure(1)
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,real(S_Sunde(:,o+(k-1)*ord)),freq,real(S_Theofilos(:,o+(k-1)*ord)))
        hold all
    end
    hold off
    subplot(ord,2,2*k)
    for o=1:ord
        semilogx(freq,imag(S_Sunde(:,o+(k-1)*ord)),freq,imag(S_Theofilos(:,o+(k-1)*ord)))
        hold all
    end
    hold off
end

figure(2)
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,abs(S_Sunde(:,o+(k-1)*ord)),freq,abs(S_Theofilos(:,o+(k-1)*ord)))
        hold all
    end
    hold off
    subplot(ord,2,2*k)
    for o=1:ord
        semilogx(freq,radtodeg(unwrap(angle(S_Sunde(:,o+(k-1)*ord)))),freq,radtodeg(unwrap(angle(S_Theofilos(:,o+(k-1)*ord)))))
        hold all
    end
    hold off
end

figure(3)
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,real(S_Ratio(:,o+(k-1)*ord)))
        hold all
    end
    hold off
    subplot(ord,2,2*k)
    for o=1:ord
        semilogx(freq,imag(S_Ratio(:,o+(k-1)*ord)))
        hold all
    end
    hold off
end

figure(4)
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,abs(S_Ratio(:,o+(k-1)*ord)))
        hold all
    end
    hold off
    subplot(ord,2,2*k)
    for o=1:ord
        semilogx(freq,unwrap(angle(S_Ratio(:,o+(k-1)*ord))))
        hold all
    end
    hold off
end