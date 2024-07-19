function [] = plotwavchar_fun(freq,ord,a,vel,Zch_mod,Ti,tag)


%% Plot attenuation and velocity
vel=abs(vel);

for o=1:ord
    lgd{o} = strcat('Mode # ',num2str(o));
end

% Calculate the required width for the legend
charWidth = 8; % approximate character width in pixels
padding = 20; % padding for the legend
maxChars = max(cellfun(@numel, lgd)); % find the longest legend entry
requiredLegendWidth = (maxChars * charWidth + padding) * ord; % total width for all legend entries

% Define initial figure width
initialWidth = 800;
figureWidth = max(initialWidth, requiredLegendWidth);

figure('Name', ['AtnConst_PhaseVel_' tag], 'Position', [100 100 figureWidth 600])

subplot(2,1,1)
for o=1:ord
    loglog(freq,a(:,o),'LineWidth',2);
    hold all
end
xlabel('Frequency [Hz]');
ylabel('\alpha [Np/m]');
title('Attenuation constant')
grid on;

% Add legend
hLegend1 = legend(lgd, 'Orientation', 'horizontal', 'Location', 'northoutside');
legendPosition1 = get(hLegend1, 'Position');
legendPosition1(1) = 0.5 - legendPosition1(3) / 2; % Center horizontally
legendPosition1(2) = .99 - legendPosition1(4); % Position just above the plot
set(hLegend1, 'Position', legendPosition1);

hold off
subplot(2,1,2)
for o=1:ord
    semilogx(freq,vel(:,o),'LineWidth',2)
    hold all
    % lgd{o} = strcat('Mode # ',num2str(o));
end
hold off
xlabel('Frequency [Hz]')
ylabel('\upsilon [m/s]')
title('Phase velocity')
% legend(lgd);
grid on

%% Plot charateristic impedance

figure('Name', ['CharImped_' tag], 'Position', [100 100 figureWidth 600])

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
grid on

% Add legend
hLegend1 = legend(lgd, 'Orientation', 'horizontal', 'Location', 'northoutside');
legendPosition1 = get(hLegend1, 'Position');
legendPosition1(1) = 0.5 - legendPosition1(3) / 2; % Center horizontally
legendPosition1(2) = .99 - legendPosition1(4); % Position just above the plot
set(hLegend1, 'Position', legendPosition1);

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
% Create a figure and tab group
fig = uifigure('Name', ['TransfMat_' tag], 'Position', [100 100 figureWidth 600]);
tg = uitabgroup(fig, 'Units', 'normalized', 'Position', [0 0 1 1] );

% Create legend labels
lgd = cell(1, ord);
for o = 1:ord
    lgd{o} = strcat('Phase # ', num2str(o));
end

% Create tabs for each mode
for k = 1:ord
    tab = uitab(tg, 'Title', ['Mode ' num2str(k)]);
    
    % Create a grid layout in the tab
    gridLayout = uigridlayout(tab, [2, 1], 'RowHeight', {'1x', '1x'}, 'ColumnWidth', {'1x'});
    
    % Create axes for magnitude plot
    ax1 = uiaxes(gridLayout);
    title(ax1, 'Transformation matrix magnitude');
    xlabel(ax1, 'Frequency [Hz]');
    ylabel(ax1, 'Magnitude');
    ax1.XScale = 'log';  % Set x-axis to logarithmic scale
    hold(ax1, 'on');
    % Add toolbar to ax1
    tb1 = axtoolbar(ax1, {'export', 'pan', 'zoomin', 'zoomout', 'restoreview'});
    
    % Create axes for angle plot
    ax2 = uiaxes(gridLayout);
    title(ax2, 'Transformation matrix angle');
    xlabel(ax2, 'Frequency [Hz]');
    ylabel(ax2, 'Angle [º]');
    ax2.XScale = 'log';  % Set x-axis to logarithmic scale
    hold(ax2, 'on');
    % Add toolbar to ax2
    tb2 = axtoolbar(ax2, {'export', 'pan', 'zoomin', 'zoomout', 'restoreview'});
    
    % Plot data
    for o = 1:ord
        semilogx(ax1, freq, abs(Ti(:, k + (o - 1) * ord)), 'LineWidth', 2);
        semilogx(ax2, freq, rad2deg(unwrap(angle(Ti(:, k + (o - 1) * ord)))), 'LineWidth', 2);
    end
    
    % Add grid
    grid(ax1, 'on');
    grid(ax2, 'on');
    box(ax1, 'on');
    box(ax2, 'on');
    
    % Add legends to both axes
    legend(ax1, lgd, 'Orientation', 'horizontal', 'Location', 'northoutside');
    % legend(ax2, lgd, 'Orientation', 'horizontal', 'Location', 'northoutside');
end


end