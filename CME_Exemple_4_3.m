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

x0 =0.05            % m
M = 10e-3           % Nm

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

figure;

[x1,t1] = impulse(Q_x*model);
[x2,t2] = impulse(Q_xM*model);

plot(t1,x1,'LineWidth',1.5); hold on;
plot(t2,x2,'LineWidth',1.5);

legend('Pretensió inicial i moment retardat','Moment amb retard','FontSize',14);
xlabel('t [s]','FontSize',14);
ylabel('x [m]','FontSize',14);

set(gca,'FontSize',14)
grid on;


%% Teorema del valor inicial i valor final
syms s_ real
F = (1)/(s_^2*m_eq+c_eq*s_+k_eq)*(M/R/s_);
limit(s_*F,s_,Inf)    % valor inicial de x(t)
limit(s_*F,s_,0)      % valor final de x(t)