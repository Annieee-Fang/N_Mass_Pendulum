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
conc = -0.25; % concavity of change of length

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
        C = .5; % amplitude of external force

        % driving frequency of external force
        omega_0 = nf+0.1;  % near natural frq
        %omega_0 = 2*nf; % far from natural frq
end

Theta_initial = [theta_initial; omega_initial];

%% ODE Solver
t_start = 0;
t_end = 15;
t_span = [t_start t_end]; 

ct0 = 0;
ct_step = 75;
ct_span = [0,ct0]; % time span of length changing period

for k = 1:ct_step
    ct = ct0 + k*(t_end-ct0)/ct_step;
    ct_span = [0,ct];

    [t, Y] = ode45(@(t,Theta) pendulum_ct(t, Theta, n, g, m, l, C, omega_0, l_max, l_min, conc, ct_span), t_span, Theta_initial);
    
    for i = 1:length(t)
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

        Yn(i,k) = Y(i,n);
        ln_func(i,k) = l(n);
    end

    t_size(k) = length(t);
end

% only keep the data before shortest t
Y_mod = zeros(min(t_size), ct_step);
L_mod = zeros(min(t_size), ct_step);
for i = 1:min(t_size)
    t_mod(i) = t(i);
    Y_mod(i,:) = Yn(i,:);
    L_mod(i,:) = ln_func(i,:);
end

N = ct_step;

%% Plot
% plotting
figure;
cmap = colormap(parula(N+1));
colors = cmap; 
hold on;
for k = 1:ct_step
    plot(t_mod, Y_mod(:,k), 'LineWidth', 1.5,'color', colors(k,:));
    title('Angle-Time Plot');
    xlabel('Time (s)');
    ylabel('Angle (radians)');
    grid on;
end

figure;
hold on;
for k = 1:ct_step
    plot(t_mod, L_mod(:,k), 'LineWidth', 1.5,'color', colors(k,:));
    title('Length-Time Plot');
    xlabel('Time (s)');
    ylabel('Length');
    grid on;
end