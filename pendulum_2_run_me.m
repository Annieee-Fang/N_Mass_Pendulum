
clear; clc; close all;

% parameters
n = 5;
g = 9.81; 
m = ones(n, 1) * 0.5;
l_0 = 0.4;
l = ones(n, 1) * l_0;
nf = sqrt(g/l_0); % natural frequency of pendulums

%% Initial Conditions
extForce_case = 'extForce';'noExtForce';
switch extForce_case
    case 'noExtForce'
        theta_initial = ones(n, 1) * (pi/6); % Initial angles
        omega_initial = zeros(n, 1); % Initial angular velocities
        C = 0; % no external force
        omega_0 = 0;
    case 'extForce'
        theta_initial = zeros(n,1); % zero initial conditions
        omega_initial = zeros(n,1);
        C = 0.5; % amplitude of external force
        omega_0 = nf+0.1 ; % driving frequency of external force 
end

Theta_initial = [theta_initial; omega_initial];

%% ODE Solver
t_start = 0;
t_end = 15;
t_span = [t_start t_end]; 

[t, Theta_solution] = ode45(@(t,Theta) pendulum_1(t, Theta, n, g, m, l, C, omega_0), t_span, Theta_initial);

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

% Plot
figure;
plot(t, theta_t, 'LineWidth', 1.5);
title('Angle-Time Plot');
xlabel('Time (s)');
ylabel('Angle (radians)');
grid on;
legend_str = cell(n, 1);
for i = 1:n
   legend_str{i} = sprintf('Theta %d', i);
end
legend(legend_str, 'Location', 'northeast');
