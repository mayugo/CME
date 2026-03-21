%% RESPOSTA TRANSITORIA MOMENT AMB RETART
% J.A.Mayugo, UdG, 2015
clear all, close all

%% Dades exercici
R = 60e-3           % m
m = 500e-3          % kg
J = 100e-3/2*R^2    % kg m2
c = 10              % Ns/m
k1 = 1              % N/m
k2 = 2              % N/m

x0 =0.075            % m
M = 5e-3            % Nm

tau = 30            % s
t_final = 50        % s

%% Sistema equivalent

m_eq=(m+J/R^2)
c_eq=c
k_eq=k1+k2

x1_0=2*k1/(k1+k2)*x0
x2_0=2*k2/(k1+k2)*x0

s = tf('s');          

Q_x0 = (- k1*x1_0 + k2*x2_0)/(s);   % moment incicial
Q_xM = M/R*exp(-30*s)/(s);          % moment retardat 30 segons
Q_x = Q_xM + Q_x0                   % moment total

%% Funció TF de Matlab en espai s

model = tf([1],[m_eq c_eq k_eq])

[x1,t1] = impulse(Q_x*model);
[x2,t2] = impulse(Q_xM*model);

figure;

plot(t1,x1,'LineWidth',1.5); hold on;
plot(t2,x2,'LineWidth',1.5);

legend('Pretensió + M(t)','només M(t)','FontSize',14);
xlabel('t [s]','FontSize',14);
ylabel('x [m]','FontSize',14);
set(gca,'FontSize',14)
grid on;

%% Teorema del valor inicial i valor final
syms s_ real

F_0   = (-k1*x1_0 + k2*x2_0)/(s_^2*m_eq + c_eq*s_ + k_eq)*(1/s_)
F_tau = (exp(-20*s_)*M/R)/(s_^2*m_eq + c_eq*s_ + k_eq)*(1/s_)

F = F_0 + F_tau;

initial_value_tau = limit(s_*F_tau,s_,Inf);     % valor inicial de x(t)
final_value_tau   = limit(s_*F_tau,s_,0);       % valor final de x(t)

initial_value = limit(s_*F,s_,Inf);    % valor inicial de x(t)
final_value   = limit(s_*F,s_,0);      % valor final de x(t)

disp('Aplicant només moment M(t)')
disp('--------------------------')
fprintf('Valor inicial sense CI: %.6f m\n', double(initial_value_tau))
fprintf('Valor final sense CI  : %.4f m\n\n', double(final_value_tau))

disp('Aplicant CI i moment M(t)')
disp('--------------------------')
fprintf('Valor inicial: %.6f m\n', double(initial_value))
fprintf('Valor final  : %.4f m\n\n', double(final_value))



