%% RESPOSTA TRANSITORIA EN CONDICIONS GENERALS NO NUL·LES
% 
clear; close all;

% -------------------------
% Dades exercici
% -------------------------
R = 60e-3;           % m
m = 500e-3;          % kg
J = 100e-3/2*R^2;    % kg·m2
c = 10;              % Ns/m
k1 = 1;              % N/m
k2 = 2;              % N/m

x0 = 0.05;           % pretensió
M = 10e-3;           % moment

t_jump  = 30;        % s
t_final = 50;        % s

% -------------------------
% Simulació
% -------------------------
[t1,x1,v1,Q] = model(m,J,R,c,k1,k2,x0,M,t_jump,t_final);
[t2,x2,v2,~] = model(m*20,J*20,R,c,k1,k2,x0,M,t_jump,t_final);

% -------------------------
% Gràfics
% -------------------------
figure;

subplot(3,1,1)
plot(t1, Q, '--')
ylabel('Q(t) [N]')
grid on
legend('Q(t) = F_0 + M(t)/r','Location','southeast')
ylim([0 0.3])

subplot(3,1,2)
plot(t1, x1, 'DisplayName','Inèrcia baixa')
hold on
plot(t2, x2, 'DisplayName','Inèrcia elevada')
ylabel('x [m]')
grid on
legend('Location','southeast')

subplot(3,1,3)
plot(t1, v1, 'DisplayName','Inèrcia baixa')
hold on
plot(t2, v2, 'DisplayName','Inèrcia elevada')
xlabel('t [s]')
ylabel('dx/dt [m/s]')
grid on

xlim([20 t_final])


function [t_full,x_full,v_full,Q_values] = model(m,J,R,c,k1,k2,x0,M,t_jump,t_final)

    % Sistema equivalent
    m_eq = m + J/R^2;
    c_eq = c;
    k_eq = k1 + k2;

    % Condicions inicials internes
    x1_0 = 2*k1/(k1+k2)*x0;
    x2_0 = 2*k2/(k1+k2)*x0;

    Q_x0 = -k1*x1_0 + k2*x2_0;

    % Força equivalent
    function val = Q(t)
        if t >= t_jump
            val = M/R + Q_x0;
        else
            val = Q_x0;
        end
    end

    % Sistema d’EDO
    function dydt = system(t,y)
        x = y(1);
        v = y(2);

        dxdt = v;
        dvdt = (Q(t) - c_eq*v - k_eq*x)/m_eq;

        dydt = [dxdt; dvdt];
    end

    % Tram 1
    tspan1 = linspace(0,t_jump,2000);
    [t1,y1] = ode45(@system, tspan1, [0;0]);

    % Tram 2
    tspan2 = linspace(t_jump,t_final,2000);
    y0_2 = y1(end,:)';
    [t2,y2] = ode45(@system, tspan2, y0_2);

    % Concatenar
    t_full = [t1; t2];
    x_full = [y1(:,1); y2(:,1)];
    v_full = [y1(:,2); y2(:,2)];

    % Q(t)
    Q_values = arrayfun(@Q, t_full);

end