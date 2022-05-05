function [] = plotZ_fun_ct(f,ord,Ztot_Carson,Ztot_Noda,Ztot_Deri,Ztot_AlDe,Ztot_Sunde,Ztot_Pettersson,Ztot_Semlyen)

%% Plot results
% Self impedance (Z11)
figure()
subplot(2,1,1)
loglog(f,squeeze(abs(Ztot_Carson(1,1,:))),f,squeeze(abs(Ztot_Noda(1,1,:))),f,squeeze(abs(Ztot_Deri(1,1,:))),f,squeeze(abs(Ztot_AlDe(1,1,:))),f,squeeze(abs(Ztot_Sunde(1,1,:))),f,squeeze(abs(Ztot_Pettersson(1,1,:))),f,squeeze(abs(Ztot_Semlyen(1,1,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
grid
title('Self impedance - Z11')

subplot(2,1,2)
loglog(f,squeeze(angle(Ztot_Carson(1,1,:))),f,squeeze(angle(Ztot_Noda(1,1,:))),f,squeeze(angle(Ztot_Deri(1,1,:))),f,squeeze(angle(Ztot_AlDe(1,1,:))),f,squeeze(angle(Ztot_Sunde(1,1,:))),f,squeeze(angle(Ztot_Pettersson(1,1,:))),f,squeeze(angle(Ztot_Semlyen(1,1,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
grid

% Mutual impedance (Z12)

number=ord;

counter=2;
figure(counter)

Zm_pet=squeeze(abs(Ztot_Pettersson(1,number,:)));

while (number-1>0)

    subplot(2,1,1)
    loglog(f,squeeze(abs(Ztot_Carson(1,number,:))),f,squeeze(abs(Ztot_Noda(1,number,:))),f,squeeze(abs(Ztot_Deri(1,number,:))),f,squeeze(abs(Ztot_AlDe(1,number,:))),f,squeeze(abs(Ztot_Sunde(1,number,:))),f,squeeze(abs(Ztot_Pettersson(1,number,:))),f,squeeze(abs(Ztot_Semlyen(1,number,:))),'LineWidth',2)
    xlabel('frequency (Hz)')
    ylabel('Magnitude (\Omega/m)')
    legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
    grid
    title(['Mutual impedance - Z1',num2str(number)])

    subplot(2,1,2)
    loglog(f,squeeze(angle(Ztot_Carson(1,number,:))),f,squeeze(angle(Ztot_Noda(1,number,:))),f,squeeze(angle(Ztot_Deri(1,number,:))),f,squeeze(angle(Ztot_AlDe(1,number,:))),f,squeeze(angle(Ztot_Sunde(1,number,:))),f,squeeze(angle(Ztot_Pettersson(1,number,:))),f,squeeze(angle(Ztot_Semlyen(1,number,:))),'LineWidth',2)
    xlabel('frequency (Hz)')
    ylabel('Angle (rad)')
    legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
    grid
    number=number-1;
    counter=counter+1;
    if counter<=ord
        figure(counter)
    end

end