%%% Obtenir l'eqüació 2 gdl en espai d'estat
% Micromechanical device de 2 gdl sotmés a una força sinusoidal
% J.A.Mayugo, UdG, 2012

clear all; close all;

m = 30e-6; % massa [kg]
c = .5e-3; % esmorteïment [Ns/m]
k = 1; % rigidesa [N/m]

x2I = 4e-6 % posició incial massa m2
v2I = 4e-3 % velocitat incial massa m2

A = [  0     0       1      0;
       0     0       0      1;
    -2*k/m  k/m   -2*c/m   c/m;
      k/m  -2*k/m   c/m  -2*c/m];

B = [ 0 ; 0; 0; 1/m];  % input una força a x_2

C = [1  0   0   0;
     0  1   0   0];    % defineix outputs x_1 i x_2

D = [0; 0];

tutorial_ss = ss(A, B, C, D)

t=[0:1e-4:0.4];

x0 = [ 0; x2I; 0; v2I];
[Y_a,T,X] = initial(tutorial_ss, x0,t); % Resposta a les condicions inicials

f=5e-6*sin(210*t);
[Y_b,T,X] = lsim(tutorial_ss,f,t);      % Resposta a la força f(t)

Z_1=Y_a(:,1)+Y_b(:,1);      % Resposta total, per superposició de les dues
Z_2=Y_a(:,2)+Y_b(:,2);

h0=figure; hold on
plot(T,f,'-','LineWidth',2);
xlabel('t [s]','FontSize',24);ylabel('Força [N]','FontSize',24)
title('Força f(t)','FontSize',24);set(gca,'FontSize',24)
saveas(h0,'Exemple_MatLab_3fplot.png')

h1=figure; hold on
plot(T,Y_a(:,1),'--','LineWidth',2);
plot(T,Y_b(:,1),'LineWidth',2);
xlabel('t [s]','FontSize',24);ylabel('posició [m]','FontSize',24)
legend('x_1(t) inicial','x_1(t) força');title('Resposta x_1(t)','FontSize',24);
set(gca,'FontSize',24)
saveas(h1,'Exemple_MatLab_3aplot.png')

h2=figure; hold on
plot(T,Y_a(:,2),'--r','LineWidth',2);
plot(T,Y_b(:,2),'r','LineWidth',2);
xlabel('t [s]','FontSize',24);ylabel('posició [m]','FontSize',24)
legend('x_2(t) inicial','x_2(t) força');title('Resposta x_2(t)','FontSize',24);
set(gca,'FontSize',24)
saveas(h2,'Exemple_MatLab_3bplot.png')

h3=figure;
hold on
plot(T,Z_1,'LineWidth',2);
plot(T,Z_2,'r','LineWidth',2);
xlabel('t [s]','FontSize',24);ylabel('posició [m]','FontSize',24)
legend('x_1(t)','x_2(t)');title('Resposta total','FontSize',24);
set(gca,'FontSize',24)
saveas(h3,'Exemple_MatLab_3cplot.png')

%figure; bode(tutorial_ss,{10^2,0.5*10^3})