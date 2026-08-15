%% Visualisation de la zone de navigation
clear; clc; close all;
init    % Lance le script init.m

% --- Figure ---
figure; hold on; axis equal; grid on;
xlim([-1, 12]); ylim([-1, 12]);
xlabel('x (m)'); ylabel('y (m)');
title('Zone de navigation — Bateau 2D');

% Obstacles (cercles pleins)
theta = linspace(0, 2*pi, 100);
for i = 1:size(obs, 1)
    cx = obs(i,1); cy = obs(i,2);
    fill(cx + R_obs*cos(theta), ...
        cy + R_obs*sin(theta), ...
        [0.8 0.2 0.2], 'EdgeColor', 'k');
end

% Départ
plot(start_pos(1), start_pos(2), 'gs', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Départ');

% Cible
plot(target_pos(1), target_pos(2), 'p', ...
    'MarkerSize', 16, 'MarkerFaceColor', 'y', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Cible');

legend('Location', 'northwest');