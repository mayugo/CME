%%% Càlcul de mecanismes - UdG
%%% J.A. Mayugo     2011
%%%
%%% Resposta temporal d'un sistema MCK
clear all; close all

m = 3.4;     % kg
c = 60;      % Ns/m
k = 27.5e3;      % N/m

x0=[0.102 0] % condicions incials [x_0, vel_0] [metres, metres/s]

%% Resposta lliure analítica
t=[0:0.001:0.7]

alpha=c/(2*m)
omg_d=sqrt(k/m-c^2/(4*m^2))

% representació 1
C_1=x0(1)
C_2=(x0(2)+alpha*x0(1))/omg_d

x1=exp(-alpha.*t).*(C_1.*cos(omg_d.*t)+C_2.*sin(omg_d.*t))

% representació 2
X=sqrt(C_1^2+C_2^2);
psi_d=atan(x0(1)*omg_d/(x0(2)+alpha*x0(1)));
Xt=X.*exp(-alpha.*t);
x2=X.*exp(-alpha.*t).*sin(omg_d.*t+psi_d);

h=figure;hold on
plot(t,x1)      % provar x1 o x2
plot(t,Xt,'--') % envolvent
xlabel('Temps [s]');ylabel('x(t) [m]'); hold off
saveas(h,'Ex_resposta_lliure.eps')

%% Resolució per espai d'estat
A = [ 0     1;
    -k/m -c/m];
B = [ 0 ; 1/m ];

C = [ 1  0 ];
D = 0;

figure;
model_ss = ss(A, B, C, D)
initial(model_ss, x0)