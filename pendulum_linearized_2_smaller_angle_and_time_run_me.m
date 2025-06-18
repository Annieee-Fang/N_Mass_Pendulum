clear; clc; close all;

% Parameters
n = 5;
g = 9.81;
m_0 = 1;
l_0 = 1;
l = ones(n, 1) * l_0;

%% Masses
mass_case = 'descent'; % Options: 'same', 'descent'
switch mass_case
    case 'same'
        m = ones(n, 1) * m_0;
    case 'descent'
        m = arrayfun(@(i) 0.1^(i-1) * m_0, 1:n)';
end

%% Initial Conditions
theta_initial = ones(n, 1) * 0.01; % smaller initial angles
omega_initial = zeros(n, 1);
Theta_initial = [theta_initial; omega_initial];

%% Build Linearized System Matrices (M and C)

M = zeros(n, n);
C = zeros(n, n);

for i = 1:n
    for j = 1:n
        mass_sum = sum(m(max(i, j):n));
        M(i, j) = mass_sum * l(i) * l(j);
    end
    C(i, i) = sum(m(i:n)) * g * l(i);
end

A_big = [zeros(n, n), eye(n, n); -inv(M) * C,  zeros(n, n)];

%% ODE Solver
t_start = 0;
t_end = 5; % changed from 15
t_span = [t_start t_end];
[t, Y] = ode45(@(t, Theta) pendulum_linearized_1(t, Theta, A_big), t_span, Theta_initial);

theta_t = Y(:, 1:n);


%-----
x = zeros(length(t), n);
y = zeros(length(t), n);
for i = 1:length(t)
    xi = 0;
    yi = 0;
    for j = 1:n
        xi = xi + l(j) * sin(theta_t(i, j));
        yi = yi - l(j) * cos(theta_t(i, j));
        x(i, j) = xi;
        y(i, j) = yi;
    end
end

x_trajectory = [zeros(length(t), 1), x];
y_trajectory = [zeros(length(t), 1), y];

%% Animation
fig = figure;
set(gcf, 'position', [476 360 560 420]);
set(gcf, 'color', 'w');

N = n;
cmap = colormap(parula(N + 1));
colors = cmap;

vidfile = VideoWriter('N-link_pendulum_linearized_animation.mp4', 'MPEG-4');
vidfile.FrameRate = 25;
vidfile.Quality = 100;
open(vidfile);

for i = 1:2:length(t)
    clf;
    hold on;
    
    for j = N:-1:1
        color_at_time = colors(j, :);
        h = plot(x_trajectory(1:i, j+1), y_trajectory(1:i, j+1), 'Color', [color_at_time, 0.2], 'LineWidth', 1.5);
    end
    
    plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    
    % strings and masses
    plot([0, x(i,1)], [0, y(i,1)], '-', 'Color', colors(1,:), 'LineWidth', 2);
    plot(x(i,1), y(i,1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'k');

    for j = 2:N
        plot([x(i,j-1), x(i,j)], [y(i,j-1), y(i,j)], '-', 'Color', colors(j,:), 'LineWidth', 2);
        plot(x(i,j), y(i,j), 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor', 'k');
    end
    
    % Figure properties
    axis equal;
    box on;
    axis_lim = sum(l);
    axis([-axis_lim, axis_lim, -axis_lim - 0.5, 0.5]);
    set(gca, 'fontsize', 14, 'ticklabelinterpreter', 'latex');
    title(sprintf('Linearized %d-link pendulum, t = %.2f s', N, t(i)), 'interpreter', 'latex', 'FontSize', 18);
    xlabel('x (m)', 'interpreter', 'latex');
    ylabel('y (m)', 'interpreter', 'latex');

    % about frames
    F = getframe(fig);
    writeVideo(vidfile, F);
end

close(vidfile);
disp('Video is completed!');

%% Angle-time plot
figure;
plot(t, theta_t, 'LineWidth', 1.5);
title('Angle-Time Plot (Linearized)', 'FontSize', 16);
xlabel('Time (s)');
ylabel('Angle (radians)');
grid on;
legend_str = cell(n, 1);
for i = 1:n
   legend_str{i} = sprintf('$\\theta_%d$', i);
end
legend(legend_str, 'Location', 'northeast', 'Interpreter', 'latex');