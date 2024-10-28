close all

f_param = figure('color','w');
f_param.Position = [1025 55 625 300];

% x = [0 5 10 15]; % deg
x = [0 2 4 6]; % mm

% % Rotation
% % y = [453.4 24.4 26.3 24.4; 453.4 341.5 345.6 310.2] / 450; % deg target pressure
% y = [227.1 451.1 445.9 400.8; 227.1 241.9 260.2 228.0] / 225; % deg ineq limit

% Translation
y = [453.4 57.9 57.8 33.5; 453.4 376.1 389.3 373.9] / 450; % deg target pressure
% y = [227.1 421.9 444.9 447.8; 227.1 263.1 265.0 247.0] / 225; % deg ineq limit

plot(x, y(1, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.2, 0.8, 0.7]);%[0.4660 0.8740 0.1880]);
hold on
plot(x, y(2, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.4, 0.5, 0.6]);
% plot(x, y(1, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.1, 0.8, 0.1]);%[0.4660 0.8740 0.1880]);
% hold on
% plot(x, y(2, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.18 0.55 0.34]);

hold off

ylim([0 1.3])
% ylim([0.95 2.45])
ylabel('Relative Deviation to Target Pressure')
% ylabel('Relative Deviation to Constrained Pressure Limit')

% xlabel('Tilting Angle in deg')
xlabel('Delta Translation in mm/sqrt(2)')

grid on

legend('Displacement without Correction', 'Displacement with Correction', 'Location', 'northeast')
% legend('Displacement without Correction', 'Displacement with Correction', 'Location', 'northwest')
fontsize(f_param, 12,"points")

% SAVE_fig_to_tikz(convertStringsToChars("lineplot"), gray());
