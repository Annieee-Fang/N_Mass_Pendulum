
clear; clc; close all;

% parameters
n = 5;
g = 9.81; 
m_0 = 1;
l_0 = 1;
l = ones(n, 1) * l_0;
l_max = 5; % range of length
l_min = 1;
nf = sqrt(g/l_0); % natural frequency of pendulums
conc = 0; % concavity of change of length

%% Masses
mass_case = 'same';'descent';
switch mass_case
    case 'same'
        m = ones(n, 1) * m_0;
    case 'descent'
        m = ones(n, 1);
        for i = 1:n
            m(i) = 0.1^(i-1) * m_0;
        end
end

%% Initial Conditions
extForce_case = 'extForce';'noExtForce';
switch extForce_case
    case 'noExtForce'
        theta_initial = ones(n, 1) * (pi/10); % Initial angles
        omega_initial = zeros(n, 1); % Initial angular velocities
        C = 0; % no external force
        omega_0 = 0;
    case 'extForce'
        theta_initial = zeros(n,1); % zero initial conditions
        omega_initial = zeros(n,1);
        C = 1; % amplitude of external force

        % driving frequency of external force
        omega_0 = nf+0.1;  % near natural frq
        %omega_0 = 2*nf; % far from natural frq
end

Theta_initial = [theta_initial; omega_initial];

%% ODE Solver
t_start = 0;
t_end = 15;
t_span = [t_start t_end]; 

ct_span = [0 10];

[t, Y] = ode45(@(t,Theta) pendulum_lt(t, Theta, n, g, m, l, C, omega_0, l_max, l_min, conc, ct_span), t_span, Theta_initial);

% Compute positions
a = zeros(length(t), n);
b = zeros(length(t), n);
for i = 1:length(t)
    xi = 0; yi = 0;
    
    % lenth
    xL = -0.02:.002:0.02;
    xR = 0.98:.002:1.02;
    order = 2;
    tol = 3e-3;
    concavity = conc;
    slope = 1;
    xL1 = -0.05:.05:0.05;

    a = [xL, 0.5-concavity+xL1, xR];
    b = [zeros(size(xL)), (0.5+concavity)+xL1*slope, ones(size(xR))];

    w = ones(size(a));
    w([1:length(xL),length(a)-length(xR)+1:length(a)]) = 10;
    sp = spaps(a,b, tol, w, order);
    %l(n) = (l_max-l_min) * fnval(sp, i/length(t)) + l_min;
    if ct_span(1)<t(i) && t(i) < ct_span(2)
        l(n) = (l_max-l_min) * fnval(sp, (t(i)-ct_span(1))/(ct_span(2)-ct_span(1))) + l_min;
    end
    if t(i)<=ct_span(1)
        l(n) = l_min;
    end
    if t(i)>=ct_span(2)
        l(n) = l_max;
    end

    for j = 1:n
        xi = xi + l(j)*sin(Y(i,j));
        yi = yi - l(j)*cos(Y(i,j));
        x(i,j) = xi;
        y(i,j) = yi;
    end

    % Save the current position for the trajectory plot
    x_trajectory(i,:) = x(i,:);
    y_trajectory(i,:) = y(i,:);
end

% include the pivot point (0,0)
x_trajectory = [zeros(size(x_trajectory,1),1), x_trajectory];
y_trajectory = [zeros(size(y_trajectory,1),1), y_trajectory];

N = n;

% plotting
fig = figure;
cmap = colormap(parula(N+1));
colors = cmap; 


% Set up the VideoWriter object
vidfile = VideoWriter('N-link_pendulum_animation.mp4', 'MPEG-4');
vidfile.FrameRate = 25;   % Slower playback if the number is smaller (was 25)
vidfile.Quality = 100;
open(vidfile);

set(gcf,'position',[476 360 400 300])
set(gcf,'color','w')
for i = 1:2:length(t)
    clf;
    hold on;

    for j = N:-1:1
        color_at_time = colors(j, :);  % fixed color for each mass
        h = plot(x_trajectory(1:i,j+1), y_trajectory(1:i,j+1), 'Color', color_at_time, 'LineWidth', 1.5);
        alpha_val = .1; % transparency for whole trajectory
        set(h, 'Color', [color_at_time, alpha_val]);  % apply transparency to the color
    end

    plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % pivot point at (0, 0)
    plot(x(i,1), y(i,1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(1,:),...
        'Color', colors(j,:),'MarkerEdgeColor','k');
    plot([0 x(i,1)],[0 y(i,1)], '-','Color', colors(1,:), 'LineWidth', 2)
    for j = 2:N
        plot([x(i,j-1), x(i,j)], [y(i,j-1), y(i,j)], '-','Color', colors(j,:), 'LineWidth', 2);
        plot(x(i,j), y(i,j), 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(j,:),...
            'Color', colors(j,:),'MarkerEdgeColor','k')
    end

    axis equal;
    box on
    set(gca,'fontsize',16,'ticklabelinterpreter','latex')
    title(sprintf('%d-link pendulum with strings, t = %.2f s', N, t(i)),'interpreter','latex','FontSize',20);
    axis([-N N -N-2 .5]);
    xticks(-6:2:6)
    yticks(-6:2:0)
    % drawnow;
    % Capture the frame and write to video
    F = getframe(fig);
    writeVideo(vidfile, F);
end

% Close the video file after the loop
close(vidfile);

disp('Video is completed!');



%% Animation
%figure;
%total_length = sum(l);
%axis_boundary = [-total_length, total_length, -total_length, 0.5];
%p_handle = plot(0,0,'o-','MarkerFaceColor','b','LineWidth',2); 
%title('N-Mass Pendulum Animation');
%xlabel('x (m)');
%ylabel('y (m)');
%grid on;
%axis equal;

%for k = 1:5:length(t)
    %
%    X_k = [0; x_t(k, :)'];
%    Y_k = [0; y_t(k, :)'];
%    set(plot(0, 0, 'o-', 'MarkerFaceColor', 'b', 'LineWidth', 2), 'XData', X_k, 'YData', Y_k);
%    axis(axis_boundary);
%    pause(0.01);
%end

% Plot
%figure;
%plot(t, theta_t, 'LineWidth', 1.5);
%title('Angle-Time Plot');
%xlabel('Time (s)');
%ylabel('Angle (radians)');
%grid on;
%legend_str = cell(n, 1);
%for i = 1:n
%   legend_str{i} = sprintf('Theta %d', i);
%end
%legend(legend_str, 'Location', 'northeast');

%% Plot
figure;
hold on;
for i = 1:n
    plot(t, Y(:,i), 'LineWidth', 1.5,'color', colors(i,:));
    title('Angle-Time Plot');
    xlabel('Time (s)');
    ylabel('Angle (radians)');
    grid on;
end
%legend(legend_str, 'Location', 'northeast');
