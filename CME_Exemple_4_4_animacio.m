function CME_jou_escoces()
% Matlab animation response of one degree of freedom system. 
% by J.A.Mayugo, UdG, 2016.

clear all; close all;

%% Dades del problema
t_start = 0;    % [s] temps inicial
t_end   = 20;   % [s] temps final 
t_eps1  = 1;    % [s] temps inicial aplicaci? for?a/moment
t_eps2  = 2;    % [s] temps final aplicaci? for?a/moment

R = 100e-3;     % [m] radi
k = 10e3;       % [N/m] rigidesa ressort  
m = 5;          % [kg] massa
J = 1/2*500*R^2;% [kg m^2] in?rcia de massa
c = 1000;       % [Ns/m] esmorte?ment

%M_o = [-42 -20 40 41 48]; % [Nm] Valor de moment aplicat
M_o = 100; % [Nm] Valor de moment aplicat

initial_position = 0;
initial_speed    = 0;

%% Runge-Kutta integration
for i=length(M_o)-0
    M_o_=M_o(i)
    [t,x]=ode23t(@rhs,[t_start t_end],[initial_position initial_speed]); 
    for ii=1:length(t)
        dx(ii,:)=rhs(t(ii),x(ii,:));
    end
end

%% Plot results
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
subplot(2,2,1);plot(t,x(:,1));
xlabel('t[s]'); ylabel('\theta [rad]');grid
bx = gca;f1 = hgtransform('Parent',bx);
hold on
plot(t(1),x(1,1),'ro','Parent',f1);
subplot(2,2,2);plot(t,x(:,2));
xlabel('t[s]'); ylabel('\omega [rad/s]');grid
cx = gca;f2 = hgtransform('Parent',cx);
hold on
plot(t(1),x(1,2),'ro','Parent',f2); % pause(1)
subplot(2,2,[3,4]);
angle=x(1,1);
xp0=R*cos(angle);yp0=R*sin(angle);
axis equal;axis([-R*2 3*R -1.2*R 1.2*R]);
ax = gca;
h1 = hgtransform('Parent',ax); % rotacio barra
h2 = hgtransform('Parent',ax); % translaci? seguidor
h3 = hgtransform('Parent',ax); % rotacio punts
hold on
plot([0;xp0],[0;yp0],'Color',[0 .6 .6],'LineWidth',10,'Parent',h1);
plot([xp0;xp0+1.5*R],[0;0],'Color',[0 0 .9],'LineWidth',10,'Parent',h2);
plot([xp0;xp0],[-R;R],'Color',[0 0 .9],'LineWidth',10,'Parent',h2);
plot(0,0,'r+','Parent',h3);
plot(xp0,yp0,'r+','Parent',h3);
TT = title(sprintf('Moment=%4.0f Nm. Time=%4.1f s, Angle=%4.0f graus',M_o_,t(1),mod(angle*180/pi,360))); 
hold off;

for i=2:length(t)
    angle=x(i,1);xp=R*cos(angle);yp=R*sin(angle);
    o1 = makehgtform('translate',t(i)-t(1),x(i,1)-x(1,1),0);
    o2 = makehgtform('translate',t(i)-t(1),x(i,2)-x(1,2),0);
    m1 = makehgtform('translate',xp-xp0,yp0,0);
    m2 = makehgtform('zrotate',angle);
    h1.Matrix = m2;
    h2.Matrix = m1;
    h3.Matrix = m2;
    f1.Matrix = o1;
    f2.Matrix = o2;
    TT.String = sprintf('Moment=%4.0f Nm. Time=%4.1f s, Angle=%4.0f graus',M_o_,t(i),mod(angle*180/pi,360));
    drawnow
end

%% RHS
function xdot=rhs(t,x) 
    % solves m x''+ m2 x'2 c x' + k x = f(t)
        %xdot_1 = x(2);
        %xdot_2 = -(m2/m) x(2)^2 -(c/m)*x(2) -(k/m)*x(1) + Q_x(t)/m; 
        %xdot = [xdot_1 ; xdot_2 ];
    m_eq = (m*R^2*sin(x(1))^2+J);  %   m_eq = (m*R^2+J)
    xdot_1 = x(2);   
    xdot_2 = -(m*R^2*sin(x(1))*cos(x(1))/m_eq)*x(2)^2 ...
             -(c*R^2*sin(x(1))^2/m_eq)*x(2)-(k*R^2*sin(x(1))*(1-cos(x(1)))/m_eq)...
             +Q_x(t)/m_eq;     
    xdot = [xdot_1 ; xdot_2 ];  % velocitat; acceleraci?
end

%% LHS
function Q=Q_x(t)       % Funci? for?a/moment aplicada/aplicat
    P = M_o_;           % force/moment amplitud
    if t<t_eps1||t>t_eps2   %impulse 
         Q=0;
    else
         Q=P;
    end
    %Q=P;                       % unit step
    %Q=P*t;                     % ramp input
    %Q=P*omega^2*sin(omega*t);  % harmonica
end
end