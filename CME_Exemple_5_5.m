clear; clc; close all;

%% Paràmetres

% Geometria
zp = 12;
zr = 36;
i = zr / zp;

R = 0.1;

% Inèrcies
J_m = 0.02;
J_p = 0.005;
J_t = 0.5;

m_c = 10.0;

% Rigideses
k_t = 160.0;
k = 40e3;

% Dissipació
c_r = 5.0;

% Càrregues
Mm = 10.0;
g = 9.8;

t_f = 1.5;

%% Matrius del sistema

A = [ ...
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
   -k_t/J_m,  i*k_t/J_m, 0, -c_r/J_m, 0, 0;
    i*k_t/(i^2*J_p + J_t), ...
   -(i^2*k_t + R^2*k)/(i^2*J_p + J_t), ...
    R*k/(i^2*J_p + J_t), ...
    0 0 0;
    0, R*k/m_c, -k/m_c, 0, 0, 0];

B = [ ...
    0 0 0;
    0 0 0;
    0 0 0;
    1/J_m 0 0;
    0 1/(i^2*J_p + J_t) 0;
    0 0 1/m_c];

C = [ ...
    0 0 1 0 0 0;
    0 0 0 0 0 1;
    0 0 0 c_r 0 0;
    0 R -1 0 0 0];

D = zeros(4,3);

sys = ss(A,B,C,D);

%% Simulació

t = linspace(0, t_f, 1000);

u = zeros(length(t),3);
u(:,1) = Mm;
u(:,3) = -m_c*g;

[y, t, x] = lsim(sys, u, t);

%% Figura

figure;

subplot(2,2,1)
plot(t, y(:,1)*1e3)
ylabel('Posicio $x(t)$ [$\times 10^{-3}$ m]', 'Interpreter', 'latex')
grid on

subplot(2,2,2)
plot(t, y(:,2))
ylabel('Velocitat $\dot{x}(t)$ [m/s]', 'Interpreter', 'latex')
grid on

subplot(2,2,3)
plot(t, y(:,3))
xlabel('Temps $t$ [s]', 'Interpreter', 'latex')
ylabel('$M_{dA}(t)$ [Nm]', 'Interpreter', 'latex' )
grid on

subplot(2,2,4)
plot(t, y(:,4)*1e3)
xlabel('Temps $t$ [s]', 'Interpreter', 'latex')
ylabel('$\Delta L(t)$ [$\times 10^{-3}$ m]', 'Interpreter', 'latex')
grid on

% mateix límit X a tots
for k_ax = 1:4
    subplot(2,2,k_ax)
    xlim([-0.02 t_f])
end

%% Règim estacionari

dot_theta0 = (Mm - m_c*g*R/i)/(c_r);
dot_x0 = dot_theta0/i*R;
M_dA0 = dot_theta0*c_r;
Delta_L = m_c*g/k;

fprintf('\n--------------------------------\n')
fprintf('VALORS EN RÈGIM ESTACIONARI\n')
fprintf('--------------------------------\n')
fprintf('Velocitat del motor: %.2f rad/s\n', dot_theta0)
fprintf('Velocitat del cable: %.3f m/s\n\n', dot_x0)
fprintf('Moment dissipat:     %.2f Nm\n', M_dA0)
fprintf('Elongació del cable: %.4f m\n', Delta_L)
fprintf('--------------------------------\n')