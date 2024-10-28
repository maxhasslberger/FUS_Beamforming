close all

f_param = figure('color','w');
f_param.Position = [1025 55 625 300];

x = [1 2 3];
y = [1 4 1; 3 6 2; 1 1 1; 3 2 3];

plot(x, y(1, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.1, 0.8, 0.1]);%[0.4660 0.8740 0.1880]);
hold on
plot(x, y(2, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.2, 0.8, 0.7]);
plot(x, y(3, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.18 0.55 0.34]);%[0.18 0.55 0.34]);
plot(x, y(4, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.4, 0.5, 0.6]);%[0.5 0.5 0.5]);
hold off

xlim([0.75 3.25])
ylim([0 7])
ylabel('# Target Points')

xticks([1 2 3])
xticklabels({'Amygdala', 'Caudate', 'Hippocampus'})
grid on

legend('TrX: 50', 'TrX: 1/sqrt2', 'TrZ: 50', 'TrZ: 1/sqrt2', 'Location', 'northwest')
fontsize(f_param, 12,"points")

% SAVE_fig_to_tikz(convertStringsToChars("lineplot"), gray());
