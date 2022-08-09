function [] = plotY_fun_ct(f,ord,Ytot_Imag,Ytot_Pet)

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

number=ord;

Ym_pet=squeeze(abs(Ytot_Pet(1,number,:)));

counter=ord+2;
figure(counter)
while (number-1>0)

        subplot(2,1,1)
        loglog(f,squeeze(abs(Ytot_Imag(1,number,:))),f,squeeze(abs(Ytot_Pet(1,number,:))),'LineWidth',2)
        xlabel('frequency (Hz)')
        ylabel('Magnitude (\Omega/m)')
        legend('Image','Pettersson')
        grid
        title(['Mutual admittance - Y1',num2str(number)])

        subplot(2,1,2)
        loglog(f,squeeze(angle(Ytot_Imag(1,number,:))),f,squeeze(angle(Ytot_Pet(1,number,:))),'LineWidth',2)
        xlabel('frequency (Hz)')
        ylabel('Angle (rad)')
        legend('Image','Pettersson')
        grid
            number=number-1;
            counter=counter+1;
                if counter<=2*ord
                     figure(counter)
                end
            
end
