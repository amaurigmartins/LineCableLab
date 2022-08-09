function [] = plotY_fun(f,Ytot_Imag,Ytot_Pet)

%% Plot results
% Self admittance (Y11)
figure()
subplot(2,1,1)
loglog(f,squeeze(abs(Ytot_Imag(1,1,:))),f,squeeze(abs(Ytot_Pet(1,1,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Image','Pettersson')
grid
title('Self admittance - Y11')

subplot(2,1,2)
loglog(f,squeeze(angle(Ytot_Imag(1,1,:))),f,squeeze(angle(Ytot_Pet(1,1,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Image','Pettersson')
grid

% Mutual admittance (Y12)
figure()
subplot(2,1,1)
loglog(f,squeeze(angle(Ytot_Imag(1,2,:))),f,squeeze(angle(Ytot_Pet(1,2,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Image','Pettersson')
grid
title('Mutual admittance - Y12')

subplot(2,1,2)
loglog(f,squeeze(angle(Ytot_Imag(1,2,:))),f,squeeze(angle(Ytot_Pet(1,2,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Image','Pettersson')
grid

% Mutual impedance (Z13)
figure()
subplot(2,1,1)
loglog(f,squeeze(abs(Ytot_Imag(1,3,:))),f,squeeze(abs(Ytot_Pet(1,3,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Magnitude (\Omega/m)')
legend('Image','Pettersson')
grid
title('Mutual admittance - Y13')

subplot(2,1,2)
loglog(f,squeeze(angle(Ytot_Imag(1,3,:))),f,squeeze(angle(Ytot_Pet(1,3,:))),'LineWidth',2)
xlabel('frequency (Hz)')
ylabel('Angle (rad)')
legend('Image','Pettersson')
grid

