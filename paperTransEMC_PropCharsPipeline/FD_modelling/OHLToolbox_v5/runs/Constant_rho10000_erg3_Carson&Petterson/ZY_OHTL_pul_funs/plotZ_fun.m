function [] = plotZ_fun(f,Ztot_Carson,Ztot_Noda,Ztot_Deri,Ztot_AlDe,Ztot_Sunde,Ztot_Pettersson,Ztot_Semlyen)

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
figure()
subplot(2,1,1)
loglog(f,squeeze(abs(Ztot_Carson(1,2,:))),f,squeeze(abs(Ztot_Noda(1,2,:))),f,squeeze(abs(Ztot_Deri(1,2,:))),f,squeeze(abs(Ztot_AlDe(1,2,:))),f,squeeze(abs(Ztot_Sunde(1,2,:))),f,squeeze(abs(Ztot_Pettersson(1,2,:))),f,squeeze(abs(Ztot_Semlyen(1,2,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
grid
title('Mutual impedance - Z12')

subplot(2,1,2)
loglog(f,squeeze(angle(Ztot_Carson(1,2,:))),f,squeeze(angle(Ztot_Noda(1,2,:))),f,squeeze(angle(Ztot_Deri(1,2,:))),f,squeeze(angle(Ztot_AlDe(1,2,:))),f,squeeze(angle(Ztot_Sunde(1,2,:))),f,squeeze(angle(Ztot_Pettersson(1,2,:))),f,squeeze(angle(Ztot_Semlyen(1,2,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
grid

% Mutual impedance (Z13)
figure()
subplot(2,1,1)
loglog(f,squeeze(abs(Ztot_Carson(1,3,:))),f,squeeze(abs(Ztot_Noda(1,3,:))),f,squeeze(abs(Ztot_Deri(1,3,:))),f,squeeze(abs(Ztot_AlDe(1,3,:))),f,squeeze(abs(Ztot_Sunde(1,3,:))),f,squeeze(abs(Ztot_Pettersson(1,3,:))),f,squeeze(abs(Ztot_Semlyen(1,3,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
grid
title('Mutual impedance - Z13')

subplot(2,1,2)
loglog(f,squeeze(angle(Ztot_Carson(1,3,:))),f,squeeze(angle(Ztot_Noda(1,3,:))),f,squeeze(angle(Ztot_Deri(1,3,:))),f,squeeze(angle(Ztot_AlDe(1,3,:))),f,squeeze(angle(Ztot_Sunde(1,3,:))),f,squeeze(angle(Ztot_Pettersson(1,3,:))),f,squeeze(angle(Ztot_Semlyen(1,3,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen')
grid

