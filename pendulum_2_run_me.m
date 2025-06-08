
clear; clc; close all;

n = 5;
g = 9.81; 
m = ones(n, 1) * 0.5;
l = ones(n, 1) * 0.4;

%% Initial Conditions

theta_initial = ones(n, 1) * (pi/6); % Initial angles

omega_initial = zeros(n, 1); % Initial angular velocities

Theta_initial = [theta_initial; omega_initial];

%% ODE Solver
t_start = 0;
t_end = 15;
t_span = [t_start t_end]; 

[t, Theta_solution] = ode45(@(t,Theta) pendulum_1(t, Theta, n, g, m, l), t_span, Theta_initial);

theta_t = Theta_solution(:, 1:n);

% Cartesian coordinates
x_t = zeros(size(theta_t, 1), n);
y_t = zeros(size(theta_t, 1), n);

for i = 1:n
    for j = 1:i
        x_t(:, i) = x_t(:, i) + l(j) * sin(theta_t(:, j));
        y_t(:, i) = y_t(:, i) - l(j) * cos(theta_t(:, j));
    end
end

% Animation
figure;
total_length = sum(l);
axis_boundary = [-total_length, total_length, -total_length, 0.5];
p_handle = plot(0,0,'o-','MarkerFaceColor','b','LineWidth',2); 
title('N-Mass Pendulum Animation');
xlabel('x (m)');
ylabel('y (m)');
grid on;
axis equal;

for k = 1:5:length(t)
    %
    X_k = [0; x_t(k, :)'];
    Y_k = [0; y_t(k, :)'];
    set(plot(0, 0, 'o-', 'MarkerFaceColor', 'b', 'LineWidth', 2), 'XData', X_k, 'YData', Y_k);
    axis(axis_boundary);
    pause(0.01);
end