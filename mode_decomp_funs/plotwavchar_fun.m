function [] = plotwavchar_fun(freq,ord,a,vel,Zch_mod,Ti,tag)



%% Plot attenuation and velocity
vel=abs(vel);
figure('Name', ['AtnConst_PhaseVel_' tag])
subplot(2,1,1)
for o=1:ord
    loglog(freq,a(:,o),'LineWidth',2);
    hold all
    lgd{o} = strcat('Mode # ',num2str(o));
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
    lgd{o} = strcat('Mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('\upsilon [m/s]')
title('Phase velocity')
% legend(lgd);
grid on

%% Plot charateristic impedance
% set_plot_params()
figure('Name', ['CharImped_' tag])

subplot(2,1,1)
for o=1:ord
    semilogx(freq,abs(Zch_mod(:,o)),'LineWidth',2)
    hold all
    lgd{o} = strcat('Mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('Magnitude [\Omega]')
title('Characteristic impedance')
legend(lgd);
grid on

subplot(2,1,2)
for o=1:ord
    semilogx(freq,rad2deg(unwrap(angle(Zch_mod(:,o)))),'LineWidth',2)
    hold all
    lgd{o} = strcat('Mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('Angle [º]')
% legend(lgd);
grid on

%% Plot modal transform matrix
%[sk col k_scaling w_corr h_corr k_width_height fnt_scaling
% set_plot_params([1 2 2 .8 .6 1 .7])
figure('Name', ['TransfMat_' tag])
%set_plot_params(.75,2.5,1.15) % a lot of trial n error to get this almost right

lgd = cell(1, ord);
for o = 1:ord
    lgd{o} = strcat('Mode # ', num2str(o));
end

for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,abs(Ti(:,k+(o-1)*ord)),'LineWidth',2)
        %semilogx(freq,abs(Ti(:,k+(o-1)*ord)),freq,abs(Ti_vf(:,k+(o-1)*ord))) % ??? vector fitting
        hold all
        % if o==1;lgd{o} = strcat('mode # ',num2str(o));end
    end
    % if k==1; legend(lgd); end
    grid on
    xlabel('Frequency [Hz]')
    if k==1;title('Transformation matrix magnitude');end
    hold off

    subplot(ord,2,2*k)
    for o=1:ord
        %semilogx(freq,rad2deg(angle(Ti(:,k+(o-1)*ord))))
        semilogx(freq,rad2deg(unwrap(angle(Ti(:,k+(o-1)*ord)))),'LineWidth',2)
        % if o==1;lgd{o} = strcat('mode # ',num2str(o));end
        %semilogx(freq,rad2deg(acos(real(Ti(:,k+(o-1)*ord))./abs(Ti(:,k+(o-1)*ord)))))
        %semilogx(freq,rad2deg(asin(imag(Ti(:,k+(o-1)*ord))./abs(Ti(:,k+(o-1)*ord)))))
        %semilogx(freq,rad2deg(unwrap(angle(Ti(:,k+(o-1)*ord)))),freq,rad2deg(unwrap(angle(Ti_vf(:,k+(o-1)*ord))))) % ??? vector fitting
        hold all
    end
    %     legend(lgd);
    grid on
    xlabel('Frequency [Hz]')
    if k==1;title('Transformation matrix angle');end
    hold off
end

% Create a new axis for the legend and place it above the subplots
hLegend = legend(lgd, 'Orientation', 'horizontal', 'Location', 'northoutside');
% Adjust the legend to the width of the plot area
newPosition = [0.1 0.95 0.8 0.05];
set(hLegend, 'Position', newPosition, 'Units', 'normalized');

% % Create figure with a specific name
%     figure('Name', ['TransfMat_' tag]);
%
%     % Predefine legend entries, assuming the user isn't a moron and ord > 0
%     lgd = cell(1, ord);
%     for o = 1:ord
%         lgd{o} = strcat('Mode # ', num2str(o));
%     end
%
%     % Loop over each order to create subplots
%     for k = 1:ord
%         % Magnitude plot
%         subplot(ord, 2, 2 * k - 1);
%         hold on; % Because we need to keep all plots on the same axis
%         for o = 1:ord
%             semilogx(freq, abs(Ti(:, k + (o - 1) * ord)), 'LineWidth', 2);
%         end
%         if k == 1
%             title('Transformation matrix magnitude');
%         end
%         grid on;
%         xlabel('Frequency [Hz]');
%         hold off;
%
%         % Angle plot
%         subplot(ord, 2, 2 * k);
%         hold on;
%         for o = 1:ord
%             semilogx(freq, rad2deg(unwrap(angle(Ti(:, k + (o - 1) * ord)))), 'LineWidth', 2);
%         end
%         if k == 1
%             title('Transformation matrix angle');
%         end
%         grid on;
%         xlabel('Frequency [Hz]');
%         hold off;
%     end
%
%     % Create a new axis for the legend and place it above the subplots
%     hLegend = legend(lgd, 'Orientation', 'horizontal', 'Location', 'northoutside');
%     % Adjust the legend to the width of the plot area
%     newPosition = [0.1 0.95 0.8 0.05];
%     set(hLegend, 'Position', newPosition, 'Units', 'normalized');


end