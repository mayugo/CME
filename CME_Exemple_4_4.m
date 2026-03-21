function CME_jou_escoces()
% Matlab animation response of one degree of freedom system. 
% by J.A.Mayugo, UdG, 2016.

clear all; close all;

%% Dades del problema
t_start = 0;    % [s] temps inicial
t_end   = 20;   % [s] temps final 
t_eps1  = 1;    % [s] temps inicial aplicacio forca/moment
t_eps2  = 2;    % [s] temps final aplicacio forca/moment

R = 100e-3;     % [m] radi
k = 10e3;       % [N/m] rigidesa ressort  
m = 5;          % [kg] massa
J = 1/2*500*R^2;% [kg m^2] in?rcia de massa
c = 1000;       % [Ns/m] esmorte?ment

M_o = 100;      % [Nm] Valor de moment aplicat

initial_position = 0;
initial_speed    = 0;

%% Runge-Kutta integration
[t,x]=ode23t(@rhs,[t_start t_end],[initial_position initial_speed]); 

%% Plot results
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.5]);
subplot(1,3,1);plot([t_start t_eps1 t_eps1 t_eps2 t_eps2 t_end ],[0 0 M_o M_o 0 0]);
xlabel('t[s]'); ylabel('M_0 [Nm]');grid
subplot(1,3,2);plot(t,x(:,1));
xlabel('t[s]'); ylabel('\theta [rad]');grid
subplot(1,3,3);plot(t,x(:,2));
xlabel('t[s]'); ylabel('\omega [rad/s]');grid

%% RHS
function xdot=rhs(t,x) 
    m_eq    = (m*R^2*sin(x(1))^2+J);  %   m_eq = (m*R^2+J)
    xdot_1  = x(2);   
    xdot_2  = -(m*R^2*sin(x(1))*cos(x(1))/m_eq)*x(2)^2 ...
              -(c*R^2*sin(x(1))^2/m_eq)*x(2)-(k*R^2*sin(x(1))*(1-cos(x(1)))/m_eq)...
              +Q_x(t)/m_eq;     
    xdot    = [xdot_1 ; xdot_2 ];  % velocitat; acceleracio
end

%% LHS
function Q=Q_x(t)      % Funcio forca/moment aplicada/aplicat
    if t<t_eps1||t>t_eps2   %impulse 
         Q=0;
    else
         Q=M_o;
    end
end
end