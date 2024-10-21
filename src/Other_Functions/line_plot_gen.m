close all

f_param = figure('color','w');
f_param.Position = [1025 55 625 300];

x = [1 2];
y = [1 4 1; 3 6 2]';

plot(x, y(1, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.4660 0.8740 0.1880]);
hold on
plot(x, y(2, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0 0]);
plot(x, y(3, :), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.18 0.55 0.34]);
hold off

xlim([0.75 2.25])
ylim([0 7])
ylabel('# Target Points')

xticks([1 2])
xticklabels({'50', '1/sqrt2'})
grid on

legend('Amygdala', 'Caudate', 'Hippocampus', 'Location', 'northwest')
fontsize(f_param, 12,"points")

% SAVE_fig_to_tikz(convertStringsToChars("lineplot"), gray());
