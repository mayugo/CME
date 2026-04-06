clear; clc; close all;

% Paràmetres
M = 1.0;
m = 0.5;
L = 0.6;

f = 1.0;
g = 9.81;

t_max = 16;
t_for = 6;

% Condicions inicials
y0 = [0; 0; pi/4; 0];  % [x, x_dot, theta, theta_dot]

% Temps
tspan = [0 t_max];

% Resolució ODE (equivalent a solve_ivp RK45)
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, y] = ode45(@(t,y) deriv(t,y,M,m,L,g,f,t_for), tspan, y0, opts);

% Variables
x = y(:,1);
x_dot = y(:,2);
theta = y(:,3);
theta_dot = y(:,4);

% Energies
T = 0.5*M.*x_dot.^2 + 0.5*m.*(x_dot.^2 + (L*theta_dot).^2 + 2*x_dot.*L.*theta_dot.*cos(theta));
V = m*g*L*(1 - cos(theta));
E = T + V;

%% FIGURA 2x2
figure('Position', [100 100 800 600])  % [x y width height] en píxels

% x(t)
subplot(2,2,1)
plot(t, x, 'b', 'LineWidth', 1.5)
ylabel('$x(t)$ [m]', 'Interpreter','latex')
grid on
xline(t_for, '--', 'LineWidth', 2)

% theta(t)
subplot(2,2,2)
plot(t, theta, 'Color', [1 0.5 0], 'LineWidth', 1.5)
xlabel('$t$ [s]', 'Interpreter','latex')
ylabel('$\theta(t)$ [rad]', 'Interpreter','latex')
grid on
xline(t_for, '--', 'LineWidth', 2)

% x_dot(t)
subplot(2,2,3)
plot(t, x_dot, 'b', 'LineWidth', 1.5)
ylabel('$\dot{x}(t)$ [m/s]', 'Interpreter','latex')
grid on
xline(t_for, '--', 'LineWidth', 2)

% theta_dot(t)
subplot(2,2,4)
plot(t, theta_dot, 'Color', [1 0.5 0], 'LineWidth', 1.5)
xlabel('$t$ [s]', 'Interpreter','latex')
ylabel('$\dot{\theta}(t)$ [rad/s]', 'Interpreter','latex')
grid on
xline(t_for, '--', 'LineWidth', 2)

%% FIGURA ENERGIES
figure;
plot(t, T, 'LineWidth', 1.5); hold on
plot(t, V, 'LineWidth', 1.5)
plot(t, E, '--', 'LineWidth', 1.5)

xline(t_for, '--', 'LineWidth', 2)

xlabel('$t$ [s]', 'Interpreter','latex')
ylabel('energia [J]', 'Interpreter','latex')
legend({'Energia cinetica $T$', 'Energia potencial $V$', 'Energia total $T+V$'}, ...
       'Interpreter','latex')

xlim([0 t_max])
grid on


function dydt = deriv(t, y, M, m, L, g, f, t_for)

    x = y(1);
    x_dot = y(2);
    theta = y(3);
    theta_dot = y(4);

    % Força
    if t <= t_for
        F = f;
    else
        F = 0.0;
    end

    cos_theta = cos(theta);
    sin_theta = sin(theta);

    % Calcul acceleracions resolent el sistema
    % ddx*(M + m) + m*L*ddtheta*cos(theta) - m*L*theta_dot^2*sin(theta) = F
    % L*ddtheta + ddx*cos(theta) + g*sin(theta) = 0
    % resolvem com a sistema lineal per ddx i ddtheta

    % Sistema lineal
    A = [M + m, m*L*cos_theta;
         cos_theta, L];

    b = [F + m*L*theta_dot^2*sin_theta;
         -g*sin_theta];

    sol = A \ b;
    ddx = sol(1);
    ddtheta = sol(2);

    dydt = [x_dot;
            ddx;
            theta_dot;
            ddtheta];
end